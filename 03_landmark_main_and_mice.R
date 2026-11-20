# Baseline-covariate MI and full-cohort / landmark regression models.

library(mice)
library(sandwich)
library(dplyr)
source("R/mi_helpers.R")

BLOCK_D <- c("age", "gender", "education", "marital", "residence", "smoking", "drinking")
BLOCK_C <- c("BMI", "HTN", "diabetes", "dyslipidemia", "heart_hist", "stroke_hist")
IMPUTABLE <- c("education", "marital", "residence", "smoking", "drinking", "BMI",
               "HTN", "diabetes", "dyslipidemia", "heart_hist", "stroke_hist")
lmk_cols <- c("ID", "communityID", "age", "BMI", "c11", "c13", "c15", "c18", "c20",
              "persistent_high_cesd", "gender", "education", "marital", "residence",
              "smoking", "drinking", "HTN", "diabetes", "dyslipidemia", "heart_hist",
              "stroke_hist", "Sleep_Label")

f_fc_m1 <- dep_high ~ Sleep_Label + age + gender
f_fc_m2 <- reformulate(c("Sleep_Label", BLOCK_D), response = "dep_high")
f_fc_m3 <- reformulate(c("Sleep_Label", BLOCK_D, BLOCK_C), response = "dep_high")
f_lmk_l0 <- persistent_high_cesd ~ Sleep_Label + age + gender
f_lmk_l1 <- reformulate(c("Sleep_Label", BLOCK_D), response = "persistent_high_cesd")
f_lmk_l2 <- reformulate(c("Sleep_Label", BLOCK_D, BLOCK_C), response = "persistent_high_cesd")
f_lmk_l3 <- update(f_lmk_l1, . ~ . + c11)
f_lmk_l4 <- update(f_lmk_l2, . ~ . + c11)
f_lmk_l5 <- update(f_lmk_l2, . ~ . + c11 + c13 + c15)

# Rubin pooling with Barnard-Rubin finite-sample degrees of freedom (dfcom = G - 1)
pool_rubin_scalar <- function(Q, U, dfcom, scale = c("identity", "exp", "pp")) {
  scale <- match.arg(scale)
  ps <- pool.scalar(Q = Q, U = U, n = dfcom + 1, k = 1, rule = "rubin1987")
  se <- sqrt(ps$t)
  t_crit <- qt(0.975, df = ps$df)
  lo <- ps$qbar - t_crit * se
  hi <- ps$qbar + t_crit * se
  p_val <- 2 * pt(-abs(ps$qbar / se), df = ps$df)

  if (scale == "exp") {
    return(tibble(
      estimate = exp(ps$qbar), ci_lower = exp(lo), ci_upper = exp(hi),
      log_est = ps$qbar, se = se, df = ps$df, dfcom = dfcom, fmi = ps$fmi, p_value = p_val
    ))
  }
  if (scale == "pp") {
    return(tibble(
      estimate = ps$qbar * 100, ci_lower = lo * 100, ci_upper = hi * 100,
      se = se * 100, df = ps$df, dfcom = dfcom, fmi = ps$fmi, p_value = p_val
    ))
  }
  tibble(
    estimate = ps$qbar, ci_lower = lo, ci_upper = hi,
    se = se, df = ps$df, dfcom = dfcom, fmi = ps$fmi, p_value = p_val
  )
}

# Within-imputation logistic regression, community-clustered HC1 variance, and
# regression standardization over the empirical covariate distribution
fit_pooled_logistic_std <- function(completed_list, formula_obj, exposure_var = "Sleep_Label",
                                    ref_level = "Longer-sleep", exp_level = "Shorter-sleep") {
  m     <- length(completed_list)
  G     <- length(unique(completed_list[[1]]$communityID))
  dfcom <- G - 1L
  coef_term <- paste0(exposure_var, exp_level)

  q_all <- u_all <- vector("list", m)
  q_logor <- u_logor <- numeric(m)
  q_p0    <- u_p0    <- numeric(m)
  q_p1    <- u_p1    <- numeric(m)
  q_rd    <- u_rd    <- numeric(m)
  q_logrr <- u_logrr <- numeric(m)

  for (j in seq_len(m)) {
    dat <- completed_list[[j]]
    stopifnot(!anyNA(dat[all.vars(formula_obj)]), !anyNA(dat$communityID))
    fit <- glm(formula_obj, data = dat, family = binomial("logit"), na.action = na.fail)
    stopifnot(isTRUE(fit$converged), nobs(fit) == nrow(dat), !anyNA(coef(fit)))
    V   <- vcovCL(fit, cluster = dat$communityID, type = "HC1", cadjust = TRUE)
    b   <- coef(fit)
    q_all[[j]] <- b
    u_all[[j]] <- diag(V)

    q_logor[j] <- unname(b[coef_term])
    u_logor[j] <- unname(V[coef_term, coef_term])

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
    g_lrr <- (g1 / p1) - (g0 / p0)

    q_p0[j]    <- p0
    u_p0[j]    <- as.numeric(t(g0) %*% V %*% g0)
    q_p1[j]    <- p1
    u_p1[j]    <- as.numeric(t(g1) %*% V %*% g1)
    q_rd[j]    <- p1 - p0
    u_rd[j]    <- as.numeric(t(g_rd) %*% V %*% g_rd)
    q_logrr[j] <- log(p1 / p0)
    u_logrr[j] <- as.numeric(t(g_lrr) %*% V %*% g_lrr)
  }

  list(
    n        = nrow(completed_list[[1]]),
    events   = sum(model.response(model.frame(formula_obj, completed_list[[1]]))),
    per_imputation = data.frame(imputation = seq_len(m), log_or = q_logor,
                                variance = u_logor, OR = exp(q_logor)),
    log_or_mcse = sqrt(var(q_logor) / m),
    coefficients = bind_rows(lapply(names(q_all[[1]]), function(term) {
      Q <- vapply(q_all, function(x) unname(x[term]), numeric(1))
      U <- vapply(u_all, function(x) unname(x[term]), numeric(1))
      mutate(pool_rubin_scalar(Q, U, dfcom, scale = "exp"), term = term, .before = 1)
    })),
    clusters = G,
    m        = m,
    or       = pool_rubin_scalar(q_logor, u_logor, dfcom, scale = "exp"),
    prob_ref = pool_rubin_scalar(q_p0,    u_p0,    dfcom, scale = "pp"),
    prob_exp = pool_rubin_scalar(q_p1,    u_p1,    dfcom, scale = "pp"),
    rd_pp    = pool_rubin_scalar(q_rd,    u_rd,    dfcom, scale = "pp"),
    rr       = pool_rubin_scalar(q_logrr, u_logrr, dfcom, scale = "exp")
  )
}


run_main_models <- function() {
  input <- "output/02_trajectories"
  out <- "output/03_main_models"
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  fc <- readRDS(file.path(input, "full_cohort_with_trajectories_16717.rds"))
  lmk <- readRDS(file.path(input, "landmark_primary_with_trajectories_7888.rds"))
  stopifnot(nrow(fc) == 16717L, nrow(lmk) == 7888L,
            sum(lmk$persistent_high_cesd) == 2051L, !anyDuplicated(fc$ID), !anyDuplicated(lmk$ID))

  fc_frame <- fc[c("ID", "communityID", BLOCK_D, BLOCK_C, "Sleep_Label", "Dep_Label")]
  rownames(fc_frame) <- NULL
  mi_fc <- run_covariate_mice(fc_frame, c(BLOCK_D, BLOCK_C, "Sleep_Label", "Dep_Label"),
                             profile = "standard")
  completed_fc <- lapply(seq_len(mi_fc$m), function(j) {
    d <- complete(mi_fc, j)
    d$dep_high <- as.integer(d$Dep_Label == "Higher-symptom")
    d
  })

  lmk_frame <- lmk[lmk_cols]
  rownames(lmk_frame) <- NULL
  mi_lmk <- run_covariate_mice(lmk_frame,
    c(BLOCK_D, BLOCK_C, "Sleep_Label", "persistent_high_cesd", "c11", "c13", "c15"), profile = "primary")
  completed_lmk <- lapply(seq_len(mi_lmk$m), function(j) complete(mi_lmk, j))
  fc_formulas <- list(Model1 = f_fc_m1, Model2 = f_fc_m2, Model3 = f_fc_m3)
  lmk_formulas <- list(L0 = f_lmk_l0, L1 = f_lmk_l1, L2 = f_lmk_l2,
                      L3 = f_lmk_l3, L4 = f_lmk_l4, L5 = f_lmk_l5)
  results <- list(
    full_cohort = lapply(fc_formulas, function(f) fit_pooled_logistic_std(completed_fc, f)),
    landmark = lapply(lmk_formulas, function(f) fit_pooled_logistic_std(completed_lmk, f)))

  saveRDS(mi_fc, file.path(out, "mids_full_cohort_16717.rds"))
  saveRDS(completed_fc, file.path(out, "completed_full_cohort_16717.rds"))
  saveRDS(mi_lmk, file.path(out, "mids_landmark_7888.rds"))
  saveRDS(completed_lmk, file.path(out, "completed_landmark_7888.rds"))
  saveRDS(results, file.path(out, "primary_model_estimates.rds"))
  write_mi_diagnostics(mi_fc, out, "full_cohort_mi")
  write_mi_diagnostics(mi_lmk, out, "landmark_mi")
  for (branch in names(results)) {
    rows <- bind_rows(lapply(names(results[[branch]]), function(model) {
      x <- results[[branch]][[model]]
      bind_rows(lapply(c("or", "prob_ref", "prob_exp", "rd_pp", "rr"), function(measure) {
        mutate(x[[measure]], model = model, measure = measure, n = x$n,
               events = x$events, clusters = x$clusters, .before = 1)
      }))
    }))
    write.csv(rows, file.path(out, paste0(branch, "_estimates.csv")), row.names = FALSE)
    coefficients <- bind_rows(lapply(names(results[[branch]]), function(model) {
      mutate(results[[branch]][[model]]$coefficients, model = model, .before = 1)
    }))
    write.csv(coefficients, file.path(out, paste0(branch, "_coefficients.csv")), row.names = FALSE)
  }
  capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
  invisible(results)
}

if (sys.nframe() == 0L) run_main_models()
