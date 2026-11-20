# Quadratic growth mixture models for full-cohort and landmark trajectories.
# Requires lcmm 2.2.2. Generated output contains participant data: keep it private.
# Launch with BLAS/OMP thread limits of 1; PSOCK workers inherit these settings.

build_long_panel <- function(df, wave_cols, wave_years, min_obs = 1L) {
  stopifnot(nrow(df) > 0L, !anyNA(df$ID), !anyDuplicated(df$ID),
            length(wave_cols) == length(wave_years), !anyDuplicated(wave_years),
            all(wave_cols %in% names(df)), all(vapply(df[wave_cols], is.numeric, logical(1))))
  df <- df[order(df$ID), , drop = FALSE]
  id_map <- data.frame(ID = df$ID, ID_num = seq_len(nrow(df)))
  long_df <- do.call(rbind, lapply(seq_along(wave_cols), function(i) {
    ok <- !is.na(df[[wave_cols[i]]])
    data.frame(ID = df$ID[ok], ID_num = as.integer(which(ok)),
               wave = rep(as.integer(wave_years[i]), sum(ok)),
               time = rep(as.numeric(wave_years[i] - 2011L), sum(ok)),
               val = as.numeric(df[[wave_cols[i]]][ok]))
  }))
  stopifnot(all(is.finite(long_df$val)), setequal(long_df$wave, wave_years),
            all(tabulate(long_df$ID_num, nrow(id_map)) >= min_obs))
  long_df <- long_df[order(long_df$ID_num, long_df$time), , drop = FALSE]
  rownames(long_df) <- rownames(id_map) <- NULL
  list(long = long_df, id_map = id_map, waves = as.integer(wave_years))
}

fit_value <- function(mod, name) {
  x <- if (inherits(mod, "hlme")) mod[[name]] else NULL
  if (is.numeric(x) && length(x) == 1L && is.finite(x)) x else NA_real_
}

finite_hlme <- function(mod, k) {
  isTRUE(fit_value(mod, "ng") == k) && is.numeric(mod$best) &&
    length(mod$best) > 0L && all(is.finite(mod$best)) && is.finite(fit_value(mod, "loglik"))
}

needs_continuation <- function(mod, k) {
  isTRUE(fit_value(mod, "conv") == 2L) && finite_hlme(mod, k)
}

trajectory_model_call <- function(k, maxiter = 500L, B = NULL) {
  args <- list(fixed = quote(val ~ time + I(time^2)), random = quote(~ time),
               subject = "ID_num", ng = k, nwg = FALSE, idiag = FALSE,
               data = as.name(".trajectory_data"), maxiter = maxiter,
               convB = 1e-4, convL = 1e-4, convG = 1e-4, nproc = 1L,
               na.action = 2L, var.time = "time", partialH = FALSE,
               verbose = FALSE, returndata = FALSE)
  if (k > 1L) args$mixture <- quote(~ time + I(time^2))
  if (!is.null(B)) args$B <- B
  as.call(c(list(as.name("hlme")), args))
}

trajectory_grid_call <- function(k) {
  substitute(lcmm::gridsearch(m = MODEL, rep = 32L, maxiter = 30L,
                              minit = minit, cl = cl),
             list(MODEL = trajectory_model_call(k)))
}

search_start <- function(k, seed, minit, cl) {
  # lcmm 2.2.2 runs the short fits on workers; intercept only its final local call.
  hlme <- function(B, ...) list(selected_B = B)
  set.seed(seed, kind = "L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  found <- eval(trajectory_grid_call(k))
  stopifnot(is.numeric(found$selected_B), length(found$selected_B) > 0L,
            all(is.finite(found$selected_B)))
  found$selected_B
}

fit_call <- function(k, maxiter = 500L, B = NULL) {
  mc <- trajectory_model_call(k, maxiter, B)
  mc[[1L]] <- quote(lcmm::hlme)
  eval(mc, envir = .GlobalEnv)
}

capture_fit <- function(expr) {
  warnings <- character()
  result <- tryCatch(list(value = withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = NULL), error = function(e) list(value = NULL, error = conditionMessage(e)))
  c(result, list(warnings = warnings))
}

posterior_quality <- function(pp, id_map, k) {
  prob_cols <- paste0("prob", seq_len(k))
  if (!is.data.frame(pp) || !setequal(names(pp), c("ID_num", "class", prob_cols)) ||
      anyDuplicated(names(pp)) || nrow(pp) != nrow(id_map) || anyDuplicated(pp$ID_num) ||
      !setequal(pp$ID_num, id_map$ID_num)) stop("Unexpected posterior layout or fitted ID set.")
  pp <- pp[match(id_map$ID_num, pp$ID_num), , drop = FALSE]
  P <- as.matrix(pp[prob_cols])
  modal <- pp$class
  if (!is.numeric(P) || any(!is.finite(P)) || any(P < 0 | P > 1) ||
      any(abs(rowSums(P) - 1) > 1e-6) || !is.numeric(modal) || anyNA(modal) ||
      any(!modal %in% seq_len(k)) || any(modal != max.col(P, ties.method = "first"))) {
    stop("Invalid posterior probabilities or modal assignments.")
  }
  counts <- tabulate(modal, k)
  avepp <- vapply(seq_len(k), function(g) {
    if (counts[g] == 0L) NA_real_ else mean(P[modal == g, g])
  }, numeric(1))
  classes <- data.frame(class = seq_len(k), n_modal = counts,
                        prop_modal = counts / nrow(id_map), AvePP = avepp)
  positive <- P[P > 0]
  entropy <- if (k == 1L) NA_real_ else 1 + sum(positive * log(positive)) / (nrow(P) * log(k))
  list(pprob = data.frame(id_map, class = as.integer(modal), pp[prob_cols], row.names = NULL),
       classes = classes, entropy = entropy,
       criteria_met = isTRUE(all(classes$prop_modal >= 0.05 & avepp >= 0.70)))
}

summarize_hlme_fit <- function(mod, id_map, k, seed = NA_integer_) {
  converged <- isTRUE(fit_value(mod, "conv") == 1L) && finite_hlme(mod, k) &&
    is.finite(fit_value(mod, "BIC"))
  row <- data.frame(k = k, seed = seed, conv = fit_value(mod, "conv"),
                    loglik = fit_value(mod, "loglik"), BIC = fit_value(mod, "BIC"),
                    converged = converged, min_prop = NA_real_, min_avepp = NA_real_,
                    entropy = NA_real_, criteria_met = FALSE, admissible = FALSE,
                    prediction_ok = FALSE, error = "")
  quality <- list(pprob = NULL, classes = NULL)
  if (converged) {
    checked <- capture_fit({
      stopifnot(fit_value(mod, "ns") == nrow(id_map))
      posterior_quality(mod$pprob, id_map, k)
    })
    if (!is.null(checked$error)) {
      row$error <- checked$error
    } else {
      quality <- checked$value
      row$min_prop <- min(quality$classes$prop_modal)
      row$min_avepp <- min(quality$classes$AvePP)
      row$entropy <- quality$entropy
      row$criteria_met <- row$admissible <- quality$criteria_met
    }
  }
  list(summary = row, pprob = quality$pprob, classes = quality$classes, curves = NULL)
}

trajectory_curves <- function(pred, k, waves) {
  cols <- if (k == 1L) "Ypred" else paste0("Ypred_class", seq_len(k))
  times <- as.numeric(waves - 2011L)
  stopifnot(is.matrix(pred$pred), is.numeric(pred$pred),
            identical(colnames(pred$pred), cols), nrow(pred$pred) == length(times),
            all(is.finite(pred$pred)), is.data.frame(pred$times),
            identical(as.numeric(pred$times$time), times))
  means <- unname(colMeans(pred$pred))
  ranks <- order(order(means, seq_len(k)))
  data.frame(class = rep(seq_len(k), each = length(times)), wave = waves, time = times,
             predicted_mean = as.vector(pred$pred),
             wave_mean = rep(means, each = length(times)),
             display_rank = rep(ranks, each = length(times)))
}

fit_candidate <- function(k, seed, minit, cl, panel_obj) {
  grid <- full <- continuation <- NULL
  if (k > 1L) grid <- capture_fit(search_start(k, seed, minit, cl))
  if (k == 1L || is.null(grid$error)) full <- capture_fit(fit_call(k, 500L, grid$value))
  if (needs_continuation(full$value, k)) {
    continuation <- capture_fit(fit_call(k, 500L, full$value$best))
  }
  mod <- if (!is.null(continuation)) continuation$value else full$value
  meta <- summarize_hlme_fit(mod, panel_obj$id_map, k, seed)
  meta$summary$full_conv <- fit_value(full$value, "conv")
  meta$summary$full_niter <- fit_value(full$value, "niter")
  meta$summary$continuation_conv <- fit_value(continuation$value, "conv")
  meta$summary$continuation_niter <- fit_value(continuation$value, "niter")
  phases <- list(grid = grid, full = full, continuation = continuation)
  errors <- c(unlist(lapply(phases, function(x) x$error)), meta$summary$error)
  meta$summary$warnings <- paste(unlist(lapply(phases, function(x) x$warnings)), collapse = " | ")
  if (meta$summary$converged) {
    curves <- capture_fit(trajectory_curves(
      lcmm::predictY(mod, newdata = data.frame(time = as.numeric(panel_obj$waves - 2011L)),
                    var.time = "time", draws = FALSE, na.action = 2L), k, panel_obj$waves))
    meta$curves <- curves$value
    meta$summary$prediction_ok <- is.null(curves$error)
    errors <- c(errors, curves$error)
    meta$summary$warnings <- paste(c(meta$summary$warnings, curves$warnings), collapse = " | ")
  }
  meta$summary$error <- paste(errors[nzchar(errors)], collapse = " | ")
  list(model = mod, grid = grid, full = full, continuation = continuation, meta = meta)
}

select_candidate <- function(status) {
  # Select the converged BIC winner within each k before testing class quality.
  best <- integer()
  for (k in sort(unique(status$k))) {
    idx <- which(status$k == k & status$converged & is.finite(status$BIC))
    if (length(idx)) best <- c(best, idx[order(status$BIC[idx], status$seed[idx], na.last = TRUE)[1L]])
  }
  eligible <- best[which(status$criteria_met[best])]
  winner <- if (!length(eligible)) integer() else
    eligible[order(status$BIC[eligible], status$k[eligible])[1L]]
  list(per_k = best, winner = winner)
}

save_grid_summary <- function(candidates, out_dir) {
  status <- do.call(rbind, lapply(names(candidates), function(id) {
    data.frame(candidate_id = id, candidates[[id]]$meta$summary, row.names = NULL)
  }))
  pick <- select_candidate(status)
  status$per_k_best <- seq_len(nrow(status)) %in% pick$per_k
  status$selected <- seq_len(nrow(status)) %in% pick$winner
  utils::write.csv(status, file.path(out_dir, "model_selection.csv"), row.names = FALSE, na = "")
  for (field in c("classes", "curves")) {
    rows <- lapply(names(candidates), function(id) {
      x <- candidates[[id]]$meta[[field]]
      if (!is.null(x)) data.frame(candidate_id = id, x, row.names = NULL)
    })
    rows <- Filter(Negate(is.null), rows)
    if (length(rows)) utils::write.csv(do.call(rbind, rows),
      file.path(out_dir, if (field == "classes") "class_quality.csv" else "fitted_curves.csv"),
      row.names = FALSE, na = "")
  }
  list(grid_table = status,
       selected = if (length(pick$winner)) candidates[[status$candidate_id[pick$winner]]] else NULL,
       candidates = candidates)
}

fit_trajectory_grid <- function(panel_obj, out_dir, max_k = 4L, seeds = 2025:2027, workers = 24L) {
  stopifnot(length(max_k) == 1L, max_k %in% 1:4,
            length(workers) == 1L, workers %in% 1:31,
            length(seeds) > 0L, !anyNA(seeds), !anyDuplicated(seeds))
  threads <- c(OMP_NUM_THREADS = "1", OMP_THREAD_LIMIT = "1", OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1", BLIS_NUM_THREADS = "1", GOTO_NUM_THREADS = "1",
               VECLIB_MAXIMUM_THREADS = "1")
  old_threads <- Sys.getenv(names(threads), unset = NA_character_)
  do.call(Sys.setenv, as.list(threads))
  on.exit({
    Sys.unsetenv(names(old_threads)[is.na(old_threads)])
    do.call(Sys.setenv, as.list(old_threads[!is.na(old_threads)]))
  }, add = TRUE)
  if (!requireNamespace("lcmm", quietly = TRUE) || utils::packageVersion("lcmm") != "2.2.2") {
    stop("Trajectory fitting requires lcmm 2.2.2.")
  }
  suppressPackageStartupMessages(library("lcmm", character.only = TRUE))
  # gridsearch's worker export resolves its data symbol in the global environment.
  if (exists(".trajectory_data", envir = .GlobalEnv, inherits = FALSE)) stop(".trajectory_data is already in use.")
  assign(".trajectory_data", panel_obj$long[c("ID_num", "time", "val")], envir = .GlobalEnv)
  on.exit(rm(list = ".trajectory_data", envir = .GlobalEnv), add = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cl <- NULL
  on.exit(if (!is.null(cl)) parallel::stopCluster(cl), add = TRUE)
  candidates <- list()
  add_candidate <- function(id, candidate) {
    candidates[[id]] <<- candidate
    saveRDS(candidate, file.path(out_dir, paste0(id, ".rds")))
    save_grid_summary(candidates, out_dir)
  }
  result <- add_candidate("k1", fit_candidate(1L, NA_integer_, NULL, NULL, panel_obj))
  minit <- candidates$k1$model
  if (!isTRUE(fit_value(minit, "conv") == 1L) || !finite_hlme(minit, 1L)) {
    stop("The one-class initializer did not converge; its fit and diagnostics were saved in ", out_dir)
  }
  if (max_k > 1L) {
    cl <- parallel::makePSOCKcluster(rep("localhost", workers), rscript_args = "--vanilla", outfile = "")
    versions <- parallel::clusterCall(cl, function(libs) {
      .libPaths(libs)
      as.character(utils::packageVersion("lcmm"))
    }, .libPaths())
    stopifnot(all(unlist(versions) == "2.2.2"))
    for (k in seq.int(2L, max_k)) for (s in seeds) {
      id <- sprintf("k%d_seed%d", k, s)
      result <- add_candidate(id, fit_candidate(k, s, minit, cl, panel_obj))
    }
  }
  result
}

assign_neutral_labels <- function(selected_fit, panel_obj, lower_label, higher_label, ref_label) {
  if (is.null(selected_fit)) stop("No eligible trajectory solution; inspect model_selection.csv.")
  if (!isTRUE(selected_fit$meta$summary$k == 2L)) {
    stop("Selected k = ", selected_fit$meta$summary$k,
         "; binary downstream analyses require k = 2. Classes have not been merged.")
  }
  curves <- selected_fit$meta$curves
  if (!isTRUE(selected_fit$meta$summary$prediction_ok) || is.null(curves)) {
    stop("Selected fitted curves are unavailable; inspect model_selection.csv.")
  }
  levels <- unique(curves[c("class", "wave_mean")])
  levels <- levels[order(levels$wave_mean, levels$class), , drop = FALSE]
  stopifnot(nrow(levels) == 2L, ref_label %in% c(lower_label, higher_label))
  pp <- selected_fit$meta$pprob
  stopifnot(identical(pp$ID, panel_obj$id_map$ID), all(pp$class %in% levels$class))
  labels <- c(lower_label, higher_label)[match(pp$class, levels$class)]
  pp$traj_label <- factor(labels, levels = c(ref_label, setdiff(c(lower_label, higher_label), ref_label)))
  pp$post_prob_assigned <- as.matrix(pp[c("prob1", "prob2")])[cbind(seq_len(nrow(pp)), pp$class)]
  pp
}

attach_trajectory_labels <- function(df, classes, label_col, prob_col = NULL) {
  idx <- match(df$ID, classes$ID)
  stopifnot(!anyNA(idx), !anyDuplicated(classes$ID))
  df[[label_col]] <- classes$traj_label[idx]
  if (!is.null(prob_col)) df[[prob_col]] <- classes$post_prob_assigned[idx]
  df
}

main <- function(in_dir = "output/01_cohort", out_dir = "output/02_trajectories",
                 workers = c(fc_sleep = 15L, fc_cesd10 = 5L, lmk_sleep = 24L,
                             paired_c10 = 14L, paired_c9 = 14L)) {
  stopifnot(all(c("fc_sleep", "fc_cesd10", "lmk_sleep", "paired_c10", "paired_c9") %in% names(workers)))
  master_df <- readRDS(file.path(in_dir, "charls_analytic_master.rds"))
  master_df <- master_df[order(master_df$ID), , drop = FALSE]
  results <- list()
  run_domain <- function(df, cols, waves, name, symptom = FALSE, min_obs = 1L, max_k = 4L) {
    panel <- build_long_panel(df, cols, waves, min_obs)
    result <- fit_trajectory_grid(panel, file.path(out_dir, name), max_k = max_k,
                                  workers = workers[[name]])
    results[[name]] <<- result$grid_table
    saveRDS(results, file.path(out_dir, "trajectory_model_selection_grids.rds"))
    assign_neutral_labels(result$selected, panel,
      lower_label = if (symptom) "Lower-symptom" else "Shorter-sleep",
      higher_label = if (symptom) "Higher-symptom" else "Longer-sleep",
      ref_label = if (symptom) "Lower-symptom" else "Longer-sleep")
  }

  waves <- c(2011L, 2013L, 2015L, 2018L, 2020L)
  sleep_cols <- c("s11", "s13", "s15", "s18", "s20")
  cesd_cols <- c("c11", "c13", "c15", "c18", "c20")
  fc_df <- master_df[which(master_df$in_full_cohort), , drop = FALSE]
  fc_sleep <- run_domain(fc_df, sleep_cols, waves, "fc_sleep")
  fc_cesd10 <- run_domain(fc_df, cesd_cols, waves, "fc_cesd10", symptom = TRUE)
  fc_analytic <- attach_trajectory_labels(fc_df, fc_sleep, "Sleep_Label", "sleep_prob_assigned")
  fc_analytic <- attach_trajectory_labels(fc_analytic, fc_cesd10, "Dep_Label", "dep_prob_assigned")
  fc_analytic$dep_high <- as.integer(fc_analytic$Dep_Label == "Higher-symptom")

  lmk_df <- master_df[which(master_df$training_eligible), , drop = FALSE]
  lmk_sleep <- run_domain(lmk_df, sleep_cols[1:3], waves[1:3], "lmk_sleep", min_obs = 2L)
  lmk_training_analytic <- attach_trajectory_labels(lmk_df, lmk_sleep, "Sleep_Label", "sleep_prob_assigned")
  lmk_primary_analytic <- lmk_training_analytic[which(lmk_training_analytic$common_eligible), , drop = FALSE]

  paired_df <- fc_analytic[which(fc_analytic$in_full_paired9), , drop = FALSE]
  paired_c10 <- run_domain(paired_df, cesd_cols, waves, "paired_c10", symptom = TRUE, max_k = 2L)
  paired_c9 <- run_domain(paired_df, c("c9_11", "c9_13", "c9_15", "c9_18", "c9_20"),
                          waves, "paired_c9", symptom = TRUE, max_k = 2L)
  fc_paired9_analytic <- attach_trajectory_labels(paired_df, paired_c10, "Dep10_Paired_Label")
  fc_paired9_analytic <- attach_trajectory_labels(fc_paired9_analytic, paired_c9, "Dep9_Label")
  fc_paired9_analytic$dep10_high <- as.integer(fc_paired9_analytic$Dep10_Paired_Label == "Higher-symptom")
  fc_paired9_analytic$dep9_high <- as.integer(fc_paired9_analytic$Dep9_Label == "Higher-symptom")

  saveRDS(fc_analytic, file.path(out_dir, "full_cohort_with_trajectories_16717.rds"))
  saveRDS(fc_paired9_analytic, file.path(out_dir, "full_cohort_paired_cesd9_16693.rds"))
  saveRDS(lmk_training_analytic, file.path(out_dir, "landmark_training_with_trajectories_13851.rds"))
  saveRDS(lmk_primary_analytic, file.path(out_dir, "landmark_primary_with_trajectories_7888.rds"))
  invisible(results)
}

if (sys.nframe() == 0L) main()
