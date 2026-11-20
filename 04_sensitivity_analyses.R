# Sensitivity analyses for sample restrictions, symptom scoring, sleep duration,
# missing covariates, survey and response weights, and baseline pain/ADL.

library(mice)
library(sandwich)
library(survey)
library(splines)
library(dplyr)

source("03_landmark_main_and_mice.R", local = TRUE)

align_completed <- function(completed, observed) {
  stopifnot(length(completed) > 1L, !anyNA(observed$ID), !anyDuplicated(observed$ID))
  lapply(completed, function(d) {
    stopifnot(!anyNA(d$ID), !anyDuplicated(d$ID), setequal(d$ID, observed$ID))
    d <- d[match(observed$ID, d$ID), , drop = FALSE]
    stopifnot(identical(d$communityID, observed$communityID))
    d
  })
}

paired_mi_frame <- function(dat) {
  dat$Dep_Label <- dat$Dep10_Paired_Label
  aux <- if (identical(dat$Dep_Label, dat$Dep9_Label)) character() else "Dep9_Label"
  out <- dat[c("ID", "communityID", BLOCK_D, BLOCK_C, "Sleep_Label", "Dep_Label", aux)]
  rownames(out) <- NULL
  out
}

weight_summary <- function(name, w) {
  stopifnot(length(w) > 0L, all(is.finite(w) & w > 0))
  q <- unname(quantile(w, c(0, .01, .05, .5, .95, .99, 1), type = 7))
  ess <- sum(w)^2 / sum(w^2)
  data.frame(weight = name, n = length(w), ESS = ess, ESS_fraction = ess / length(w),
             mean = mean(w), min = q[1], p01 = q[2], p05 = q[3], p50 = q[4],
             p95 = q[5], p99 = q[6], max = q[7], max_over_mean = q[7] / mean(w),
             top_1pct_weight_share = sum(sort(w, decreasing = TRUE)[seq_len(ceiling(.01 * length(w)))]) / sum(w))
}

response_balance <- function(dat, formula_obj, raw, truncated) {
  x <- model.matrix(formula_obj, dat)[, -1L, drop = FALSE]
  observed <- x[dat$responded, , drop = FALSE]
  stopifnot(nrow(x) == nrow(dat), nrow(observed) == length(raw), length(raw) == length(truncated))
  target <- colMeans(x)
  scale <- apply(x, 2L, sd)
  scale[scale == 0] <- NA_real_
  u <- colMeans(observed)
  a <- as.numeric(crossprod(raw, observed) / sum(raw))
  b <- as.numeric(crossprod(truncated, observed) / sum(truncated))
  data.frame(term = colnames(x), target_mean = target, unweighted_response_mean = u,
             raw_response_mean = a, truncated_response_mean = b, target_sd = scale,
             unweighted_smd = (u - target) / scale, raw_smd = (a - target) / scale,
             truncated_smd = (b - target) / scale, row.names = NULL)
}

# Flatten aggregate results only; participant-level completed datasets are not exported.
estimate_rows <- function(x, analysis = "", metadata = list()) {
  if (is.data.frame(x)) {
    for (nm in setdiff(names(metadata), names(x))) x[[nm]] <- metadata[[nm]]
    return(cbind(analysis = analysis, x))
  }
  if (!is.list(x)) return(data.frame(analysis = analysis, estimate = x))
  for (nm in intersect(c("n", "events", "clusters", "m"), names(x))) metadata[[nm]] <- x[[nm]]
  children <- setdiff(names(x), c("n", "events", "clusters", "m"))
  bind_rows(lapply(children, function(nm) {
    estimate_rows(x[[nm]], paste(c(analysis[nzchar(analysis)], nm), collapse = "/"), metadata)
  }))
}

# Single-dataset logistic regression with community-clustered HC1 variance and standardization
fit_single_logistic_std <- function(dat, formula_obj, exposure_var = "Sleep_Label",
                                    ref_level = "Longer-sleep", exp_level = "Shorter-sleep") {
  G     <- length(unique(dat$communityID))
  df_cl <- G - 1L
  t_crit <- qt(0.975, df = df_cl)
  coef_term <- paste0(exposure_var, exp_level)

  stopifnot(G > 1L, !anyNA(dat$communityID))
  fit <- glm(formula_obj, data = dat, family = binomial("logit"), na.action = na.fail)
  V   <- vcovCL(fit, cluster = dat$communityID, type = "HC1", cadjust = TRUE)
  b   <- coef(fit)

  est_logor <- unname(b[coef_term])
  se_logor  <- sqrt(unname(V[coef_term, coef_term]))

  d0 <- d1 <- dat
  d0[[exposure_var]] <- factor(ref_level, levels = c(ref_level, exp_level))
  d1[[exposure_var]] <- factor(exp_level, levels = c(ref_level, exp_level))
  X0 <- model.matrix(delete.response(terms(fit)), data = d0)
  X1 <- model.matrix(delete.response(terms(fit)), data = d1)

  p0_i <- as.numeric(plogis(X0 %*% b))
  p1_i <- as.numeric(plogis(X1 %*% b))
  p0   <- mean(p0_i)
  p1   <- mean(p1_i)

  g0    <- colMeans((p0_i * (1 - p0_i)) * X0)
  g1    <- colMeans((p1_i * (1 - p1_i)) * X1)
  g_rd  <- g1 - g0
  se_rd <- sqrt(as.numeric(t(g_rd) %*% V %*% g_rd))

  tibble(
    n        = nrow(dat),
    events   = sum(model.response(model.frame(fit))),
    clusters = G,
    OR       = exp(est_logor),
    OR_lower = exp(est_logor - t_crit * se_logor),
    OR_upper = exp(est_logor + t_crit * se_logor),
    p_value  = 2 * pt(-abs(est_logor / se_logor), df = df_cl),
    RD_pp    = (p1 - p0) * 100,
    RD_lower = ((p1 - p0) - t_crit * se_rd) * 100,
    RD_upper = ((p1 - p0) + t_crit * se_rd) * 100
  )
}

# Joint stacked cluster-level score influence function comparison (Pipper et al. 2012)
fit_paired_history_comparison <- function(comp_subset, f_c9, f_c10 = f_lmk_l5) {
  m     <- length(comp_subset)
  G     <- length(unique(comp_subset[[1]]$communityID))
  dfcom <- G - 1L
  coef_term <- "Sleep_LabelShorter-sleep"

  q_dlogor <- u_dlogor <- numeric(m)
  q_drd    <- u_drd    <- numeric(m)

  for (j in seq_len(m)) {
    dat   <- comp_subset[[j]]
    n     <- nrow(dat)
    fit10 <- glm(f_c10, data = dat, family = binomial("logit"), na.action = na.fail)
    fit9  <- glm(f_c9,  data = dat, family = binomial("logit"), na.action = na.fail)

    p10 <- length(coef(fit10))
    p9  <- length(coef(fit9))
    c10_adj <- sqrt((G / (G - 1)) * ((n - 1) / (n - p10)))
    c9_adj  <- sqrt((G / (G - 1)) * ((n - 1) / (n - p9)))

    psi10 <- (estfun(fit10) %*% (bread(fit10) / n)) * c10_adj
    psi9  <- (estfun(fit9)  %*% (bread(fit9)  / n)) * c9_adj
    cl_f  <- factor(dat$communityID)
    U10   <- rowsum(psi10, group = cl_f, reorder = FALSE)
    U9    <- rowsum(psi9,  group = cl_f, reorder = FALSE)
    V_jnt <- crossprod(cbind(U10, U9))

    idx10 <- which(names(coef(fit10)) == coef_term)
    idx9  <- p10 + which(names(coef(fit9)) == coef_term)
    c_vec <- numeric(p10 + p9)
    c_vec[idx10] <- -1
    c_vec[idx9]  <-  1

    q_dlogor[j] <- unname(coef(fit9)[coef_term] - coef(fit10)[coef_term])
    u_dlogor[j] <- as.numeric(t(c_vec) %*% V_jnt %*% c_vec)

    std_grad <- function(fit_obj) {
      b  <- coef(fit_obj)
      d0 <- d1 <- dat
      d0$Sleep_Label <- factor("Longer-sleep",  levels = c("Longer-sleep", "Shorter-sleep"))
      d1$Sleep_Label <- factor("Shorter-sleep", levels = c("Longer-sleep", "Shorter-sleep"))
      X0 <- model.matrix(delete.response(terms(fit_obj)), data = d0)
      X1 <- model.matrix(delete.response(terms(fit_obj)), data = d1)
      p0_i <- as.numeric(plogis(X0 %*% b))
      p1_i <- as.numeric(plogis(X1 %*% b))
      g0   <- colMeans((p0_i * (1 - p0_i)) * X0)
      g1   <- colMeans((p1_i * (1 - p1_i)) * X1)
      list(rd = mean(p1_i) - mean(p0_i), grad = g1 - g0)
    }

    s10 <- std_grad(fit10)
    s9  <- std_grad(fit9)
    g_drd <- c(-s10$grad, s9$grad)

    q_drd[j] <- s9$rd - s10$rd
    u_drd[j] <- as.numeric(t(g_drd) %*% V_jnt %*% g_drd)
  }

  list(
    n           = nrow(comp_subset[[1]]),
    clusters    = G,
    fit_hist10  = fit_pooled_logistic_std(comp_subset, f_c10),
    fit_hist9   = fit_pooled_logistic_std(comp_subset, f_c9),
    delta_logor = pool_rubin_scalar(q_dlogor, u_dlogor, dfcom, scale = "identity"),
    or_ratio    = pool_rubin_scalar(q_dlogor, u_dlogor, dfcom, scale = "exp"),
    delta_rd_pp = pool_rubin_scalar(q_drd,    u_drd,    dfcom, scale = "pp")
  )
}

# Continuous 2018-2020 follow-up mean depressive symptom score linear regressions
fit_pooled_linear <- function(comp_subset, formula_obj, coef_term = "Sleep_LabelShorter-sleep") {
  m     <- length(comp_subset)
  G     <- length(unique(comp_subset[[1]]$communityID))
  dfcom <- G - 1L
  q_b   <- u_b <- numeric(m)
  for (j in seq_len(m)) {
    dat   <- comp_subset[[j]]
    fit_j <- lm(formula_obj, data = dat, na.action = na.fail)
    V_j   <- vcovCL(fit_j, cluster = dat$communityID, type = "HC1", cadjust = TRUE)
    q_b[j] <- unname(coef(fit_j)[coef_term])
    u_b[j] <- unname(V_j[coef_term, coef_term])
  }
  mutate(pool_rubin_scalar(q_b, u_b, dfcom, scale = "identity"),
         n = nrow(comp_subset[[1]]), clusters = G, m = m, .before = 1)
}

run_sensitivity_analyses <- function(IN_COHORT = "output/01_cohort",
                                     IN_TRAJ = "output/02_trajectories",
                                     IN_MAIN = "output/03_main_models",
                                     OUT_DIR = "output/04_sensitivity") {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  fc_df         <- readRDS(file.path(IN_TRAJ, "full_cohort_with_trajectories_16717.rds"))
  fc_paired9_df <- readRDS(file.path(IN_TRAJ, "full_cohort_paired_cesd9_16693.rds"))
  lmk_train_df  <- readRDS(file.path(IN_TRAJ, "landmark_training_with_trajectories_13851.rds"))
  lmk_df        <- readRDS(file.path(IN_TRAJ, "landmark_primary_with_trajectories_7888.rds"))
  comp_fc       <- readRDS(file.path(IN_MAIN, "completed_full_cohort_16717.rds"))
  comp_lmk      <- readRDS(file.path(IN_MAIN, "completed_landmark_7888.rds"))
  comp_fc <- align_completed(comp_fc, fc_df)
  comp_lmk <- align_completed(comp_lmk, lmk_df)

  # Full-cohort restrictions and paired nine-item trajectories
  idx_fc_base_lt10 <- which(!is.na(fc_df$c11) & fc_df$c11 < 10)
  idx_fc_sleep2w   <- which(fc_df$n_sleep_5w >= 2L)
  idx_fc_sleep3w   <- which(fc_df$n_sleep_5w >= 3L)

  fit_fc_m3_base_lt10 <- fit_pooled_logistic_std(lapply(comp_fc, `[`, idx_fc_base_lt10, ), f_fc_m3)
  fit_fc_m3_sleep2w   <- fit_pooled_logistic_std(lapply(comp_fc, `[`, idx_fc_sleep2w,   ), f_fc_m3)
  fit_fc_m3_sleep3w   <- fit_pooled_logistic_std(lapply(comp_fc, `[`, idx_fc_sleep3w,   ), f_fc_m3)

  # Paired 10-item vs 9-item symptom trajectory MI branch (n = 16,693)
  fc_paired_mi_frame <- paired_mi_frame(fc_paired9_df)

  mids_fc_paired9 <- run_covariate_mice(
    fc_paired_mi_frame,
    active_predictors = setdiff(names(fc_paired_mi_frame), c("ID", "communityID")),
    m = 100L, maxit = 20L, seed = 2025L, profile = "standard"
  )
  write_mi_diagnostics(mids_fc_paired9, OUT_DIR, "paired_full_cohort_mi")
  comp_fc_paired9 <- lapply(seq_len(100L), function(j) {
    d_j <- complete(mids_fc_paired9, action = j)
    d_j$Dep9_Label <- fc_paired9_df$Dep9_Label
    d_j$dep10_high <- as.integer(d_j$Dep_Label  == "Higher-symptom")
    d_j$dep9_high  <- as.integer(d_j$Dep9_Label == "Higher-symptom")
    d_j
  })

  f_fc_m3_cesd10 <- update(f_fc_m3, dep10_high ~ .)
  f_fc_m3_cesd9  <- update(f_fc_m3, dep9_high  ~ .)

  fit_fc_m3_cesd10_paired <- fit_pooled_logistic_std(comp_fc_paired9, f_fc_m3_cesd10)
  fit_fc_m3_cesd9_paired  <- fit_pooled_logistic_std(comp_fc_paired9, f_fc_m3_cesd9)

  # Classification agreement between 10-item and 9-item symptom trajectories (Cohen's kappa & ARI)
  tab_10v9 <- table(fc_paired9_df$Dep10_Paired_Label, fc_paired9_df$Dep9_Label)
  n_pair   <- sum(tab_10v9)
  p_obs    <- sum(diag(tab_10v9)) / n_pair
  p_exp    <- sum(rowSums(tab_10v9) * colSums(tab_10v9)) / (n_pair^2)
  kappa_10v9 <- (p_obs - p_exp) / (1 - p_exp)

  comb2  <- function(x) ifelse(x >= 2, x * (x - 1) / 2, 0)
  sum_ij <- sum(comb2(tab_10v9))
  sum_i  <- sum(comb2(rowSums(tab_10v9)))
  sum_j  <- sum(comb2(colSums(tab_10v9)))
  e_idx  <- sum_i * sum_j / comb2(n_pair)
  ari_10v9 <- (sum_ij - e_idx) / (0.5 * (sum_i + sum_j) - e_idx)

  # Landmark complete cases (n = 4,404) and below-threshold subgroup (n = 3,457)
  lmk_mi_frame <- lmk_df[lmk_cols]
  lmk_cc_df <- lmk_mi_frame[complete.cases(lmk_mi_frame), , drop = FALSE]

  fit_lmk_cc_l2 <- fit_single_logistic_std(lmk_cc_df, f_lmk_l2)
  fit_lmk_cc_l4 <- fit_single_logistic_std(lmk_cc_df, f_lmk_l4)
  fit_lmk_cc_l5 <- fit_single_logistic_std(lmk_cc_df, f_lmk_l5)

  idx_sub10 <- which(lmk_df$c11 < 10 & lmk_df$c13 < 10 & lmk_df$c15 < 10)
  comp_lmk_sub10 <- lapply(comp_lmk, `[`, idx_sub10, )
  fit_lmk_sub10_l5_mi <- fit_pooled_logistic_std(comp_lmk_sub10, f_lmk_l5)
  fit_lmk_sub10_l5_cc <- fit_single_logistic_std(
    filter(lmk_cc_df, c11 < 10 & c13 < 10 & c15 < 10), f_lmk_l5
  )

  # Continuous 2011-2015 mean sleep: linear and natural cubic spline
  knots_4   <- as.numeric(quantile(lmk_df$mean_sleep_11_15, probs = c(0.05, 0.35, 0.65, 0.95), type = 7))
  rng_99    <- as.numeric(quantile(lmk_df$mean_sleep_11_15, probs = c(0.01, 0.99), type = 7))
  ns_basis  <- ns(lmk_df$mean_sleep_11_15, knots = knots_4[2:3], Boundary.knots = knots_4[c(1, 4)], intercept = FALSE)
  colnames(ns_basis) <- paste0("sleep_ns", 1:3)

  cont_sleep_df <- cbind(lmk_df, as.data.frame(ns_basis))
  cont_mi_frame <- cont_sleep_df[c(BLOCK_D, BLOCK_C, "Sleep_Label", "c11", "c13", "c15",
                                   "persistent_high_cesd", "sleep_ns1", "sleep_ns2", "sleep_ns3")]
  rownames(cont_mi_frame) <- NULL

  mids_cont_sleep <- run_covariate_mice(
    cont_mi_frame,
    active_predictors = names(cont_mi_frame),
    m = 100L, maxit = 20L, seed = 2025L, profile = "sensitivity"
  )
  write_mi_diagnostics(mids_cont_sleep, OUT_DIR, "continuous_sleep_mi")
  comp_cont_sleep <- lapply(seq_len(100L), function(j) {
    d_j <- complete(mids_cont_sleep, action = j)
    d_j$communityID      <- lmk_df$communityID
    d_j$mean_sleep_11_15 <- lmk_df$mean_sleep_11_15
    d_j
  })

  f_lmk_sleep_lin <- persistent_high_cesd ~ mean_sleep_11_15 + age + gender +
    education + marital + residence + smoking + drinking +
    BMI + HTN + diabetes + dyslipidemia + heart_hist + stroke_hist + c11 + c13 + c15

  dfcom_lmk <- length(unique(lmk_df$communityID)) - 1L
  q_lin <- u_lin <- numeric(100L)
  for (j in seq_len(100L)) {
    fit_j <- glm(f_lmk_sleep_lin, data = comp_cont_sleep[[j]], family = binomial("logit"))
    V_j   <- vcovCL(fit_j, cluster = comp_cont_sleep[[j]]$communityID, type = "HC1", cadjust = TRUE)
    q_lin[j] <- unname(coef(fit_j)["mean_sleep_11_15"])
    u_lin[j] <- unname(V_j["mean_sleep_11_15", "mean_sleep_11_15"])
  }
  fit_lmk_sleep_linear <- pool_rubin_scalar(q_lin, u_lin, dfcom_lmk, scale = "exp")

  eval_spline_contrasts <- function(target_hours, ref_hour = 7.0) {
    b_ref  <- as.numeric(predict(ns_basis, newx = ref_hour))
    B_targ <- predict(ns_basis, newx = target_hours)
    D_mat  <- sweep(B_targ, 2L, b_ref, FUN = "-")
    colnames(D_mat) <- paste0("sleep_ns", 1:3)

    Q_mat <- U_mat <- matrix(NA_real_, nrow = length(comp_cont_sleep), ncol = length(target_hours))
    f_sp  <- persistent_high_cesd ~ sleep_ns1 + sleep_ns2 + sleep_ns3 +
      age + gender + education + marital + residence + smoking + drinking +
      BMI + HTN + diabetes + dyslipidemia + heart_hist + stroke_hist + c11 + c13 + c15

    for (j in seq_along(comp_cont_sleep)) {
      dat   <- comp_cont_sleep[[j]]
      fit_j <- glm(f_sp, data = dat, family = binomial("logit"))
      V_j   <- vcovCL(fit_j, cluster = dat$communityID, type = "HC1", cadjust = TRUE)
      sp_cols <- paste0("sleep_ns", 1:3)
      b_sp  <- coef(fit_j)[sp_cols]
      V_sp  <- V_j[sp_cols, sp_cols, drop = FALSE]

      Q_mat[j, ] <- as.numeric(D_mat %*% b_sp)
      U_mat[j, ] <- rowSums((D_mat %*% V_sp) * D_mat)
    }

    bind_rows(lapply(seq_along(target_hours), function(k) {
      result <- if (all(D_mat[k, ] == 0)) {
        tibble(estimate = 1, ci_lower = 1, ci_upper = 1, log_est = 0, se = 0,
               df = NA_real_, dfcom = dfcom_lmk, fmi = NA_real_, p_value = NA_real_)
      } else pool_rubin_scalar(Q_mat[, k], U_mat[, k], dfcom_lmk, scale = "exp")
      mutate(result, mean_sleep_h = target_hours[k], reference_h = ref_hour, .before = 1)
    }))
  }

  fit_lmk_sleep_spline <- list(
    contrasts = eval_spline_contrasts(c(5.0, 6.0, 8.0), ref_hour = 7.0),
    curve     = eval_spline_contrasts(sort(unique(c(seq(rng_99[1], rng_99[2], length.out = 161L),
                                                    5, 6, 7, 8))), ref_hour = 7.0)
  )

  # Landmark sleep-item deletion and paired history comparison
  nine_sub <- lmk_df %>%
    filter(paired8_eligible) %>%
    mutate(
      all9_complete = as.integer(paired9_eligible),
      mean9         = mean9_1820,
      mean10        = mean10_1820
    )

  nine_mi_frame <- nine_sub[c(BLOCK_D, BLOCK_C, "Sleep_Label", "c11", "c13", "c15",
                              "persistent_high_cesd", "c9_11", "c9_13", "c9_15",
                              "mean9", "mean10", "all9_complete")]
  rownames(nine_mi_frame) <- NULL

  mids_nine <- run_covariate_mice(
    nine_mi_frame,
    active_predictors = names(nine_mi_frame),
    m = 100L, maxit = 20L, seed = 2025L, profile = "sensitivity"
  )
  write_mi_diagnostics(mids_nine, OUT_DIR, "nine_item_mi")
  comp_lmk_8of9 <- lapply(seq_len(100L), function(j) {
    d_j <- complete(mids_nine, action = j)
    d_j$ID          <- nine_sub$ID
    d_j$communityID <- nine_sub$communityID
    d_j$c9strict_11 <- nine_sub$c9strict_11
    d_j$c9strict_13 <- nine_sub$c9strict_13
    d_j$c9strict_15 <- nine_sub$c9strict_15
    d_j$mean9_strict <- nine_sub$mean9_strict_1820
    d_j
  })
  comp_lmk_9of9 <- lapply(comp_lmk_8of9, function(d_j) d_j[d_j$all9_complete == 1L, ])

  f_lmk_hist3w_c10 <- persistent_high_cesd ~ Sleep_Label + age + gender +
    education + marital + residence + smoking + drinking +
    BMI + HTN + diabetes + dyslipidemia + heart_hist + stroke_hist + c11 + c13 + c15

  f_lmk_hist3w_c9 <- persistent_high_cesd ~ Sleep_Label + age + gender +
    education + marital + residence + smoking + drinking +
    BMI + HTN + diabetes + dyslipidemia + heart_hist + stroke_hist + c9_11 + c9_13 + c9_15

  fit_lmk_m3_hist3w <- fit_pooled_logistic_std(comp_lmk_8of9, f_lmk_hist3w_c10)
  fit_lmk_m3_cesd9  <- fit_pooled_logistic_std(comp_lmk_8of9, f_lmk_hist3w_c9)

  f_lmk_hist3w_c9strict <- update(f_lmk_hist3w_c9, . ~ . - c9_11 - c9_13 - c9_15 + c9strict_11 + c9strict_13 + c9strict_15)
  paired_history_8of9   <- fit_paired_history_comparison(comp_lmk_8of9, f_c9 = f_lmk_hist3w_c9)
  paired_history_9of9   <- fit_paired_history_comparison(comp_lmk_9of9, f_c9 = f_lmk_hist3w_c9strict)

  f_cont10_hist10 <- update(f_lmk_hist3w_c10, mean10 ~ .)
  f_cont10_hist9  <- update(f_lmk_hist3w_c9,  mean10 ~ .)
  f_cont9_hist9   <- update(f_lmk_hist3w_c9,  mean9  ~ .)

  fit_lmk_cont10_hist10 <- fit_pooled_linear(comp_lmk_8of9, f_cont10_hist10)
  fit_lmk_cont10_hist9  <- fit_pooled_linear(comp_lmk_8of9, f_cont10_hist9)
  fit_lmk_cont9_hist9   <- fit_pooled_linear(comp_lmk_8of9, f_cont9_hist9)
  # Table S19: history coverage on the same nine-item outcome and completions.
  f_cont9_hist2011 <- update(f_cont9_hist9, . ~ . - c9_13 - c9_15)
  fit_lmk_cont9_hist2011 <- fit_pooled_linear(comp_lmk_8of9, f_cont9_hist2011)
  history_coverage <- bind_rows(
    mutate(fit_lmk_cont9_hist2011, model = "A_2011_history", .before = 1),
    mutate(fit_lmk_cont9_hist9, model = "B_three_wave_history", .before = 1)
  )
  write.csv(history_coverage, file.path(OUT_DIR, "Table_S19_history_coverage.csv"),
            row.names = FALSE)
  f_cont9_hist9strict <- update(f_lmk_hist3w_c9strict, mean9_strict ~ .)
  fit_lmk_cont10_hist10_9of9 <- fit_pooled_linear(comp_lmk_9of9, f_cont10_hist10)
  fit_lmk_cont9_hist9_9of9 <- fit_pooled_linear(comp_lmk_9of9, f_cont9_hist9strict)

  # Baseline survey and survivor-response weighting
  # Construct the design before subsetting to the landmark domain.
  all_wt <- readRDS(file.path(IN_COHORT, "baseline_survey_design.rds")) %>%
    mutate(in_lmk_domain = ID %in% lmk_df$ID)
  all_wt <- as.data.frame(all_wt)
  stopifnot(!anyNA(all_wt$ID), !anyDuplicated(all_wt$ID), !anyNA(all_wt$communityID),
            all(is.finite(all_wt$base_weight) & all_wt$base_weight > 0))
  rownames(all_wt) <- all_wt$ID

  des_all    <- svydesign(ids = ~communityID, weights = ~base_weight, data = all_wt)
  des_domain <- subset(des_all, in_lmk_domain)
  df_svy     <- degf(des_domain)

  q_svy <- u_svy <- numeric(length(comp_lmk))
  for (j in seq_along(comp_lmk)) {
    z_j <- comp_lmk[[j]]
    p_idx <- match(des_domain$variables$ID, z_j$ID)
    stopifnot(!anyNA(p_idx), identical(des_domain$variables$communityID, z_j$communityID[p_idx]))
    des_j <- des_domain
    des_j$variables[all.vars(f_lmk_l5)] <- z_j[p_idx, all.vars(f_lmk_l5), drop = FALSE]
    fit_s <- svyglm(f_lmk_l5, design = des_j, family = quasibinomial(), na.action = na.fail)
    q_svy[j] <- unname(coef(fit_s)["Sleep_LabelShorter-sleep"])
    u_svy[j] <- unname(vcov(fit_s)["Sleep_LabelShorter-sleep", "Sleep_LabelShorter-sleep"])
  }
  fit_lmk_l5_svy <- pool_rubin_scalar(q_svy, u_svy, df_svy, scale = "exp")

  # Survivor-response IPW targeting confirmed survivors with complete 2011-2015 histories (n = 9,297)
  surv_target <- lmk_train_df %>%
    filter(survivor_resp_target) %>%
    mutate(responded = outcome_complete)

  f_resp   <- responded ~ Sleep_Label + age + gender + c11 + c13 + c15
  fit_resp <- glm(f_resp, data = surv_target, family = binomial("logit"), na.action = na.fail)
  p_resp   <- unname(fitted(fit_resp))
  stopifnot(!anyNA(surv_target$responded), all(is.finite(p_resp) & p_resp > 0 & p_resp < 1))
  sw_raw   <- mean(surv_target$responded) / p_resp[surv_target$responded]
  sw_cuts  <- unname(quantile(sw_raw, probs = c(0.01, 0.99), type = 7))
  sw_trunc <- pmax(sw_cuts[1], pmin(sw_raw, sw_cuts[2]))
  resp_ids <- surv_target$ID[surv_target$responded]

  fit_weighted_l5_series <- function(wt_vec, ids = resp_ids) {
    q_w <- u_w <- numeric(length(comp_lmk))
    for (j in seq_along(comp_lmk)) {
      idx <- match(ids, comp_lmk[[j]]$ID)
      stopifnot(!anyNA(idx), length(wt_vec) == length(ids))
      z_j <- comp_lmk[[j]][idx, , drop = FALSE]
      z_j$.w <- wt_vec
      fit_w <- glm(f_lmk_l5, data = z_j, weights = .w, family = quasibinomial(), na.action = na.fail)
      V_w   <- vcovCL(fit_w, cluster = z_j$communityID, type = "HC1", cadjust = TRUE)
      q_w[j] <- unname(coef(fit_w)["Sleep_LabelShorter-sleep"])
      u_w[j] <- unname(V_w["Sleep_LabelShorter-sleep", "Sleep_LabelShorter-sleep"])
    }
    pool_rubin_scalar(q_w, u_w, length(unique(z_j$communityID)) - 1L, scale = "exp")
  }

  fit_lmk_l5_svy_matched <- fit_weighted_l5_series(rep(1, nrow(des_domain$variables)),
                                                  ids = des_domain$variables$ID)
  fit_lmk_l5_ipw_matched <- fit_weighted_l5_series(rep(1, length(resp_ids)))
  fit_lmk_l5_ipw_raw   <- fit_weighted_l5_series(sw_raw)
  fit_lmk_l5_ipw_trunc <- fit_weighted_l5_series(sw_trunc)
  weight_diagnostics <- bind_rows(
    weight_summary("Survey matched unweighted", rep(1, nrow(des_domain$variables))),
    weight_summary("Baseline survey weight", des_domain$variables$base_weight),
    weight_summary("Response matched unweighted", rep(1, length(resp_ids))),
    weight_summary("Response IPW", sw_raw),
    weight_summary("Response IPW winsorized 1/99", sw_trunc)
  )
  balance <- response_balance(surv_target, f_resp, sw_raw, sw_trunc)
  response_diagnostics <- data.frame(
    target_n = nrow(surv_target), respondents = length(resp_ids),
    response_rate = mean(surv_target$responded), p_min = min(p_resp), p_max = max(p_resp),
    winsor_p01 = sw_cuts[1], winsor_p99 = sw_cuts[2],
    n_clipped_low = sum(sw_raw < sw_cuts[1]), n_clipped_high = sum(sw_raw > sw_cuts[2])
  )

  # Broader-health complete cases: pain and basic ADL (n = 2,847)
  lmk_health_cc <- lmk_df %>%
    filter(complete.cases(across(all_of(c(BLOCK_D, BLOCK_C, "c11", "c13", "c15", "pain", "adl_limit")))))

  fit_lmk_health_l2 <- list(
    base     = fit_single_logistic_std(lmk_health_cc, f_lmk_l2),
    expanded = fit_single_logistic_std(lmk_health_cc, update(f_lmk_l2, . ~ . + pain + adl_limit))
  )
  fit_lmk_health_l4 <- list(
    base     = fit_single_logistic_std(lmk_health_cc, f_lmk_l4),
    expanded = fit_single_logistic_std(lmk_health_cc, update(f_lmk_l4, . ~ . + pain + adl_limit))
  )
  fit_lmk_health_l5 <- list(
    base     = fit_single_logistic_std(lmk_health_cc, f_lmk_l5),
    expanded = fit_single_logistic_std(lmk_health_cc, update(f_lmk_l5, . ~ . + pain + adl_limit))
  )

  estimates <- list(
      fc_restrictions = list(
        base_lt10 = fit_fc_m3_base_lt10,
        sleep_2w  = fit_fc_m3_sleep2w,
        sleep_3w  = fit_fc_m3_sleep3w
      ),
      fc_paired_cesd9 = list(
        cesd10 = fit_fc_m3_cesd10_paired,
        cesd9  = fit_fc_m3_cesd9_paired,
        kappa  = kappa_10v9,
        ari    = ari_10v9
      ),
      lmk_complete_case = list(L2 = fit_lmk_cc_l2, L4 = fit_lmk_cc_l4, L5 = fit_lmk_cc_l5),
      lmk_sub10         = list(MI = fit_lmk_sub10_l5_mi, CC = fit_lmk_sub10_l5_cc),
      lmk_cont_sleep    = list(linear = fit_lmk_sleep_linear, spline = fit_lmk_sleep_spline),
      lmk_item_deletion = list(
        fit_lmk_m3_hist3w = fit_lmk_m3_hist3w,
        fit_lmk_m3_cesd9  = fit_lmk_m3_cesd9,
        paired_8of9       = paired_history_8of9,
        paired_9of9       = paired_history_9of9,
        cont10_hist10     = fit_lmk_cont10_hist10,
        cont10_hist9      = fit_lmk_cont10_hist9,
        cont9_hist9       = fit_lmk_cont9_hist9,
        cont10_hist10_9of9 = fit_lmk_cont10_hist10_9of9,
        cont9_hist9_9of9   = fit_lmk_cont9_hist9_9of9
      ),
      lmk_weights = list(
        survey_matched_l5 = fit_lmk_l5_svy_matched,
        survey_l5    = fit_lmk_l5_svy,
        ipw_matched_l5 = fit_lmk_l5_ipw_matched,
        ipw_raw_l5   = fit_lmk_l5_ipw_raw,
        ipw_trunc_l5 = fit_lmk_l5_ipw_trunc
      ),
      lmk_broader_health = list(
        L2 = fit_lmk_health_l2,
        L4 = fit_lmk_health_l4,
        L5 = fit_lmk_health_l5
      )
  )
  saveRDS(estimates, file.path(OUT_DIR, "sensitivity_analysis_estimates.rds"))
  write.csv(estimate_rows(estimates), file.path(OUT_DIR, "sensitivity_analysis_estimates.csv"),
            row.names = FALSE, na = "")
  write.csv(weight_diagnostics, file.path(OUT_DIR, "Table_S11_weight_diagnostics.csv"), row.names = FALSE)
  write.csv(balance, file.path(OUT_DIR, "Table_S11_response_balance.csv"), row.names = FALSE, na = "")
  write.csv(response_diagnostics, file.path(OUT_DIR, "Table_S11_response_diagnostics.csv"), row.names = FALSE)
  binary <- bind_rows(
    mutate(paired_history_8of9$fit_hist10$or, sample = "8-of-9", history = "10-item"),
    mutate(paired_history_8of9$fit_hist9$or, sample = "8-of-9", history = "9-item"),
    mutate(paired_history_9of9$fit_hist10$or, sample = "9-of-9", history = "10-item"),
    mutate(paired_history_9of9$fit_hist9$or, sample = "9-of-9", history = "9-item")
  )
  binary$n <- rep(c(nrow(nine_sub), nrow(comp_lmk_9of9[[1]])), each = 2L)
  binary$events <- rep(c(sum(nine_sub$persistent_high_cesd),
                         sum(comp_lmk_9of9[[1]]$persistent_high_cesd)), each = 2L)
  continuous <- bind_rows(
    mutate(fit_lmk_cont10_hist10, sample = "8-of-9", outcome = "mean10", history = "10-item", score_max = 30),
    mutate(fit_lmk_cont9_hist9, sample = "8-of-9", outcome = "mean9", history = "9-item", score_max = 27),
    mutate(fit_lmk_cont10_hist10_9of9, sample = "9-of-9", outcome = "mean10", history = "10-item", score_max = 30),
    mutate(fit_lmk_cont9_hist9_9of9, sample = "9-of-9", outcome = "mean9_strict", history = "9-item", score_max = 27)
  ) %>% mutate(estimate_pct_of_max = 100 * estimate / score_max,
               lower_pct_of_max = 100 * ci_lower / score_max,
               upper_pct_of_max = 100 * ci_upper / score_max)
  write.csv(binary, file.path(OUT_DIR, "Table_S15_binary.csv"), row.names = FALSE)
  write.csv(continuous, file.path(OUT_DIR, "Table_S15_continuous.csv"), row.names = FALSE)
  write.csv(estimate_rows(list(paired_8of9 = paired_history_8of9, paired_9of9 = paired_history_9of9)),
            file.path(OUT_DIR, "paired_history_comparisons.csv"), row.names = FALSE, na = "")
  invisible(estimates)
}

if (sys.nframe() == 0L) run_sensitivity_analyses()
