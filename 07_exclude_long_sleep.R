# Table S18: post hoc exclusion of mean observed 2011-2015 sleep >8 hours.
# Reuse the primary classes and all 100 saved imputations; never re-impute here.
source("03_landmark_main_and_mice.R", local = TRUE)
source("R/cluster_variance.R", local = TRUE)
source("R/standardize_risk.R", local = TRUE)

long_sleep_restriction <- function(d, observed) {
  stopifnot(!anyNA(d$ID), !anyDuplicated(d$ID),
            !anyNA(observed$ID), !anyDuplicated(observed$ID))
  i <- match(d$ID, observed$ID)
  stopifnot(!anyNA(i))
  sleep <- observed[i, c("s11", "s13", "s15")]
  values <- as.matrix(sleep)[!is.na(sleep)]
  stopifnot(all(rowSums(!is.na(sleep)) >= 2L), all(is.finite(values)),
            all(values >= 2 & values <= 12))
  mean_sleep <- rowMeans(sleep, na.rm = TRUE)
  list(mean_sleep = mean_sleep, keep = mean_sleep <= 8)
}

check_saved_completion <- function(a, d) {
  stopifnot(identical(names(a), names(d)), identical(a$ID, d$ID), !anyNA(a))
  for (v in names(d)) {
    observed <- !is.na(d[[v]])
    stopifnot(identical(a[[v]][observed], d[[v]][observed]))
  }
  invisible(TRUE)
}

long_sleep_flow <- function(d, restriction) {
  pair <- levels(d$Sleep_Label)
  scope <- list(Original = rep(TRUE, nrow(d)), Retained_mean_le8 = restriction$keep,
                Excluded_mean_gt8 = !restriction$keep)
  do.call(rbind, lapply(names(scope), function(s) {
    do.call(rbind, lapply(c("All", pair), function(g) {
      i <- scope[[s]] & (g == "All" | as.character(d$Sleep_Label) == g)
      data.frame(sample = s, sleep_class = g, n = sum(i),
        events = sum(d$persistent_high_cesd[i]),
        observed_percent = if (sum(i)) mean(d$persistent_high_cesd[i]) * 100 else NA_real_,
        communities = length(unique(d$communityID[i])),
        mean_sleep_h = if (sum(i)) mean(restriction$mean_sleep[i]) else NA_real_,
        exactly_8h = sum(restriction$mean_sleep[i] == 8))
    }))
  }))
}

fit_long_sleep_model <- function(a, formula, pair) {
  rownames(a) <- a$ID
  warnings <- character()
  fit <- withCallingHandlers(
    glm(formula, data = a, family = binomial("logit"), na.action = na.fail),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  stopifnot(isTRUE(fit$converged), !fit$boundary, !anyNA(coef(fit)), !length(warnings),
            nobs(fit) == nrow(a), fit$rank == length(coef(fit)),
            identical(rownames(model.frame(fit)), a$ID))
  z <- charls_standardize_risk(fit, a, "Sleep_Label", pair, setNames(a$communityID, a$ID))
  term <- paste0("Sleep_Label", pair[2])
  q <- c(log_or = unname(coef(fit)[term]), z$estimates)
  u <- c(log_or = unname(z$vcov_beta[term, term]), setNames(diag(z$vcov), names(z$estimates)))
  stopifnot(all(is.finite(q)), all(is.finite(u) & u >= 0))
  data.frame(n = nrow(a), events = sum(a$persistent_high_cesd),
    clusters = length(unique(a$communityID)), term = names(q), estimate = q,
    within_variance = u, warnings = length(warnings), iterations = fit$iter, row.names = NULL)
}

pool_long_sleep_results <- function(per) {
  keys <- unique(per[c("sample", "model_id", "term")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    z <- per[per$sample == k$sample & per$model_id == k$model_id & per$term == k$term, ]
    stopifnot(nrow(z) == 100L, identical(sort(z$imputation), 1:100),
              length(unique(z$clusters)) == 1L)
    t <- charls_rubin_scalar(z$estimate, z$within_variance, k$term,
                            dfcom = z$clusters[1] - 1L)$table
    analysis_estimate <- t$estimate
    if (k$term %in% c("log_or", "log_rr")) {
      t[c("estimate", "lower", "upper")] <- lapply(t[c("estimate", "lower", "upper")], exp)
    } else {
      t[c("estimate", "lower", "upper")] <- lapply(t[c("estimate", "lower", "upper")], function(x) 100 * x)
    }
    data.frame(sample = k$sample, model_id = k$model_id, n = z$n[1], events = z$events[1],
      communities = z$clusters[1], imputations = 100L, t,
      analysis_estimate = analysis_estimate, mcse_analysis = sqrt(var(z$estimate) / 100), row.names = NULL)
  }))
}

run_long_sleep_sensitivity <- function(mode = c("check", "run"),
    mi_path = "output/03_main_models/mids_landmark_7888.rds",
    observed_path = "output/02_trajectories/landmark_primary_with_trajectories_7888.rds",
    out = "output/07_long_sleep") {
  mode <- match.arg(mode)
  stopifnot(getRversion() == "4.5.1", as.character(packageVersion("mice")) == "3.19.0",
            packageVersion("sandwich") == package_version("3.1-3"))
  if (!all(file.exists(c(mi_path, observed_path))))
    stop("Saved primary MI and observed sleep inputs are required. No re-imputation is performed.")
  old_options <- options(contrasts = c("contr.treatment", "contr.poly"), na.action = "na.fail")
  on.exit(options(old_options), add = TRUE)
  mi <- readRDS(mi_path)
  observed <- readRDS(observed_path)
  if (!is.data.frame(observed)) observed <- observed$observed_frame
  stopifnot(inherits(mi, "mids"), is.data.frame(observed))
  d <- mi$data
  pair <- levels(d$Sleep_Label)
  # Accept both legacy and standardized sleep-trajectory labels.
  stopifnot(identical(pair, c("Normal Sleep", "Short Sleep")) ||
            identical(pair, c("Longer-sleep", "Shorter-sleep")))
  stopifnot(identical(names(d), lmk_cols), nrow(d) == 7888L,
            sum(d$persistent_high_cesd) == 2051L,
            mi$m == 100L, mi$iteration == 20L, mi$seed == 2025L,
            as.character(mi$version) == "3.19.0", mi$lastSeedValue[1] == 10403L,
            identical(unname(mi$where), unname(is.na(d))), NROW(mi$loggedEvents) == 0L)
  protected <- setdiff(names(d), IMPUTABLE)
  stopifnot(!anyNA(d[protected]), all(mi$method[protected] == ""),
            all(d$persistent_high_cesd == as.integer(d$c18 >= 10 & d$c20 >= 10)))
  restriction <- long_sleep_restriction(d, observed)
  i <- match(d$ID, observed$ID)
  for (v in intersect(protected, names(observed))) {
    stopifnot(identical(d[[v]], observed[[v]][i]))
  }
  keep <- restriction$keep
  stopifnot(sum(keep) == 7259L, sum(d$persistent_high_cesd[keep]) == 1929L,
            sum(!keep) == 629L, sum(!keep & d$Sleep_Label == pair[2]) == 0L,
            sum(restriction$mean_sleep == 8) == 387L,
            all(keep[restriction$mean_sleep == 8]))
  for (j in seq_len(mi$m)) check_saved_completion(mice::complete(mi, j), d)
  spec <- data.frame(model_id = c("M02_L0", "M02_L1", "M02_L2", "M03_L3", "M03_L4", "M02_L5"),
                     stringsAsFactors = FALSE)
  formulas <- list(f_lmk_l0, f_lmk_l1, f_lmk_l2, f_lmk_l3, f_lmk_l4, f_lmk_l5)
  flow <- long_sleep_flow(d, restriction)
  print(flow, row.names = FALSE)
  if (mode == "check") return(invisible(flow))

  # Explicit run only: 600 restricted fits and 100 original-L5 control fits.
  # Serial execution is portable and draws no new random numbers.
  per <- do.call(rbind, lapply(1:100, function(j) {
    a <- mice::complete(mi, j)
    r <- lapply(1:6, function(k) {
      cbind(sample = "Retained_mean_le8", model_id = spec$model_id[k], imputation = j,
            fit_long_sleep_model(a[keep, , drop = FALSE], formulas[[k]], pair))
    })
    r[[7]] <- cbind(sample = "Original_control", model_id = spec$model_id[6], imputation = j,
                    fit_long_sleep_model(a, formulas[[6]], pair))
    do.call(rbind, r)
  }))
  pooled <- pool_long_sleep_results(per)
  control <- pooled[pooled$sample == "Original_control", ]
  l5_ref <- data.frame(
    term = c("log_or", "risk_ref", "risk_exp", "rd", "log_rr"),
    estimate_original = c(1.07591654728391, 25.7728974366978, 26.7721644989423, 0.99926706224448, 1.03877183025089),
    lower_original    = c(0.923207504235072, 24.6381350917841, 24.8222672013234, -1.11004638452383, 0.959419469036428),
    upper_original    = c(1.25388540648665, 26.9076597816115, 28.7220617965612, 3.10858050901279, 1.12468732410287),
    stringsAsFactors = FALSE
  )
  old <- l5_ref[match(control$term, l5_ref$term), ]
  error <- as.matrix(control[c("estimate", "lower", "upper")]) -
    as.matrix(old[c("estimate_original", "lower_original", "upper_original")])
  if (anyNA(error) || max(abs(error)) >= 1e-7)
    stop("Original L5 does not reproduce the published reference; restricted results not written.")
  restricted <- pooled[pooled$sample == "Retained_mean_le8", ]
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  write.csv(flow, file.path(out, "long_sleep_cohort_flow.csv"), row.names = FALSE)
  write.csv(restricted, file.path(out, "long_sleep_restricted_results.csv"), row.names = FALSE)
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  invisible(restricted)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (!length(args)) args <- "check"
  stopifnot(length(args) %in% c(1L, 3L))
  if (length(args) == 1L) run_long_sleep_sensitivity(args[1]) else
    run_long_sleep_sensitivity(args[1], mi_path = args[2], observed_path = args[3])
}
