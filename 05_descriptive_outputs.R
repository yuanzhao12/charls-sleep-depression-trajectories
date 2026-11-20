# Observed-data summaries for Tables 1, S3, S4, S7, S10, S13 and S17.
# Run from the repository root after scripts 01-04. No models are refitted here.

source("04_sensitivity_analyses.R", local = TRUE)

class_correspondence <- function(dat) {
  stopifnot(!anyNA(dat[c("Sleep_Label", "Dep_Label", "communityID")]),
            all(dat$Sleep_Label %in% c("Longer-sleep", "Shorter-sleep")),
            all(dat$Dep_Label %in% c("Lower-symptom", "Higher-symptom")))
  s <- as.integer(dat$Sleep_Label == "Shorter-sleep")
  y <- as.integer(dat$Dep_Label == "Higher-symptom")
  G <- length(unique(dat$communityID))
  stopifnot(G > 1L)
  one <- function(name, num, den) {
    pr <- if (sum(den) > 0) sum(num) / sum(den) else NA_real_
    variance <- NA_real_
    ci <- c(NA_real_, NA_real_)
    if (is.finite(pr)) {
      # Ratio-estimator influence functions on the full sample, summed by community.
      influence <- (num - pr * den) / sum(den)
      u <- rowsum(influence, dat$communityID, reorder = FALSE)
      variance <- sum(u^2) * G / (G - 1)
      if (pr > 0 && pr < 1) {
        se_logit <- sqrt(variance) / (pr * (1 - pr))
        ci <- plogis(qlogis(pr) + c(-1, 1) * qt(.975, G - 1) * se_logit)
      }
    }
    data.frame(metric = name, numerator = sum(num), denominator = sum(den),
               estimate = pr, lower = ci[1], upper = ci[2],
               se_probability = sqrt(variance), communities = G, df = G - 1L)
  }
  rbind(one("Sensitivity", s * y, y), one("Specificity", (1 - s) * (1 - y), 1 - y),
        one("PPV", s * y, s), one("NPV", (1 - s) * (1 - y), 1 - s))
}

observed_characteristics <- function(dat, variables, groups, comparison = NULL) {
  bind_rows(lapply(variables, function(v) {
    x <- dat[[v]]
    lev <- if (is.numeric(x)) "Mean (SD)" else if (is.factor(x)) levels(x) else sort(unique(na.omit(x)))
    smd <- rep(NA_real_, length(lev))
    if (!is.null(comparison)) {
      a <- x[groups[[comparison[1]]]]
      b <- x[groups[[comparison[2]]]]
      if (is.numeric(x)) {
        scale <- sqrt((var(a, na.rm = TRUE) + var(b, na.rm = TRUE)) / 2)
        smd <- if (is.finite(scale) && scale > 0) (mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE)) / scale else NA_real_
      } else {
        pa <- vapply(lev, function(z) mean(a == z, na.rm = TRUE), numeric(1))
        pb <- vapply(lev, function(z) mean(b == z, na.rm = TRUE), numeric(1))
        scale <- sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
        smd <- ifelse(is.finite(scale) & scale > 0, (pa - pb) / scale, NA_real_)
      }
    }
    bind_rows(lapply(names(groups), function(nm) {
      z <- x[groups[[nm]]]
      obs <- z[!is.na(z)]
      numeric_var <- is.numeric(x)
      counts <- if (numeric_var) NA_integer_ else as.integer(table(factor(obs, levels = lev)))
      data.frame(variable = v, level = lev, group = nm, n = length(z), observed = length(obs),
                 missing = sum(is.na(z)), n_level = counts,
                 percent_observed = if (!numeric_var && length(obs)) 100 * counts / length(obs) else NA_real_,
                 mean = if (numeric_var && length(obs)) mean(obs) else NA_real_,
                 sd = if (numeric_var && length(obs) > 1L) sd(obs) else NA_real_,
                 smd = smd)
    }))
  }))
}

wave_counts <- function(counts, classes, domain) {
  groups <- c(list(Overall = seq_along(counts)), split(seq_along(counts), classes))
  bind_rows(lapply(names(groups), function(nm) {
    z <- counts[groups[[nm]]]
    n <- as.integer(table(factor(z, levels = 1:5)))
    data.frame(domain = domain, group = nm, waves = 1:5, n = n,
               denominator = length(z), percent = 100 * n / length(z))
  }))
}

sleep_distribution <- function(dat, columns, sample) {
  means <- rowMeans(dat[columns], na.rm = TRUE)
  bind_rows(lapply(levels(dat$Sleep_Label), function(nm) {
    x <- means[dat$Sleep_Label == nm & is.finite(means)]
    q <- unname(quantile(x, c(0, .25, .5, .75, 1), type = 7))
    data.frame(sample = sample, class = nm, n = length(x), mean = mean(x), sd = sd(x),
               min = q[1], q25 = q[2], median = q[3], q75 = q[4], max = q[5],
               n_lt6 = sum(x < 6), n_6to8 = sum(x >= 6 & x <= 8), n_gt8 = sum(x > 8),
               n_lt7 = sum(x < 7), percent_lt6 = 100 * mean(x < 6),
               percent_6to8 = 100 * mean(x >= 6 & x <= 8), percent_gt8 = 100 * mean(x > 8),
               percent_lt7 = 100 * mean(x < 7))
  }))
}

posterior_quality <- function(classes, assigned_probability, sample) {
  p <- assigned_probability
  stopifnot(nlevels(classes) == 2L, length(classes) == length(p), !anyNA(classes),
            all(is.finite(p) & p >= 0 & p <= 1))
  xlogx <- function(x) ifelse(x > 0, x * log(x), 0)
  entropy <- 1 + mean(xlogx(p) + xlogx(1 - p)) / log(2)
  bind_rows(lapply(levels(classes), function(nm) {
    z <- p[classes == nm]
    data.frame(sample = sample, class = nm, n = length(z), percent = 100 * length(z) / length(p),
               AvePP = mean(z), median_posterior = median(z), n_lt70 = sum(z < .7),
               n_lt80 = sum(z < .8), percent_lt70 = 100 * mean(z < .7),
               percent_lt80 = 100 * mean(z < .8), relative_entropy = entropy)
  }))
}

missingness_summary <- function(dat, sample) {
  variables <- c(BLOCK_D, BLOCK_C, "c11", "c13", "c15", "c18", "c20", "Sleep_Label")
  nmiss <- c(sum(!complete.cases(dat[c(BLOCK_D, BLOCK_C)])), colSums(is.na(dat[variables])))
  data.frame(sample = sample, variable = c("Any baseline model covariate", variables),
             n = nrow(dat), missing = unname(nmiss), percent_missing = 100 * unname(nmiss) / nrow(dat))
}

mi_result_summary <- function(result, mi, completed, analysis) {
  draws <- result$per_imputation
  if (is.null(draws)) stop("Saved model results do not contain per-imputation estimates.")
  m <- result$m
  stopifnot(is.data.frame(draws), all(c("imputation", "log_or", "variance", "OR") %in% names(draws)),
            m > 1L, nrow(draws) == m, !anyDuplicated(draws$imputation),
            setequal(draws$imputation, seq_len(m)), mi$m == m, length(completed) == m,
            all(is.finite(draws$log_or)), all(is.finite(draws$variance) & draws$variance >= 0),
            isTRUE(all.equal(draws$OR, exp(draws$log_or), check.attributes = FALSE)),
            nrow(result$or) == 1L, all(c("fmi", "estimate", "se", "df", "dfcom") %in% names(result$or)),
            nrow(result$coefficients) > 0L, "fmi" %in% names(result$coefficients))
  mcse <- sqrt(var(draws$log_or) / m)
  stopifnot(isTRUE(all.equal(result$log_or_mcse, mcse)))
  observed <- lapply(mi$data, function(x) which(!is.na(x)))
  preserved <- vapply(completed, function(x) {
    stopifnot(nrow(x) == nrow(mi$data), all(names(mi$data) %in% names(x)))
    all(vapply(names(mi$data), function(v) {
      identical(x[[v]][observed[[v]]], mi$data[[v]][observed[[v]]])
    }, logical(1)))
  }, logical(1))
  data.frame(
    analysis = analysis, completed_datasets = m, iterations = mi$iteration, seed = mi$seed,
    covariate_targets = sum(nzchar(mi$method)), recorded_mice_events = NROW(mi$loggedEvents),
    residual_missing_cells = sum(vapply(completed, function(x) sum(is.na(x)), numeric(1))),
    observed_protected_preserved = all(preserved), sleep_coefficient_fmi = result$or$fmi,
    largest_coefficient_fmi = max(result$coefficients$fmi),
    per_completion_OR_min = min(draws$OR), per_completion_OR_max = max(draws$OR),
    pooled_OR = result$or$estimate, pooled_log_OR_se = result$or$se,
    df = result$or$df, dfcom = result$or$dfcom, log_OR_mcse = mcse,
    within_imputation_variance = mean(draws$variance), between_imputation_variance = var(draws$log_or)
  )
}

cohort_flow <- function(d) {
  known <- !is.na(d$age) & !is.na(d$gender)
  ageok <- known & d$age >= 45
  train <- d$training_eligible
  groups <- list(
    linked_baseline = rep(TRUE, nrow(d)), missing_age_or_sex = !known,
    age_below45 = known & d$age < 45, age45_valid_demographics = ageok,
    no_sleep_any_wave = ageok & d$n_sleep_5w == 0,
    sleep_but_no_cesd10 = ageok & d$n_sleep_5w > 0 & d$n_cesd10_5w == 0,
    full_cohort = d$in_full_cohort, paired_full_cohort = d$in_full_paired9,
    age_sleep_eligible = d$age_sleep_eligible,
    death_by_index = d$age_sleep_eligible & d$vital_index == "DEAD_BY_INDEX",
    index_survival_unknown = d$age_sleep_eligible & d$vital_index == "UNKNOWN",
    index_conflict = d$age_sleep_eligible & d$vital_index == "CONFLICT",
    landmark_training = train, training_without_pair = train & !d$outcome_complete,
    training_with_pair = train & d$outcome_complete,
    outcome_pair_incomplete_history = train & d$outcome_complete & !d$history_complete,
    landmark_primary = d$common_eligible, paired_8of9 = d$paired8_eligible,
    paired_9of9 = d$paired9_eligible
  )
  bind_rows(lapply(names(groups), function(nm) {
    z <- groups[[nm]]
    stopifnot(!anyNA(z))
    data.frame(stage = nm, n = sum(z), events = if (nm %in% c("landmark_primary", "paired_8of9", "paired_9of9"))
      sum(d$persistent_high_cesd[z]) else NA_integer_)
  }))
}

followup_counts <- function(d) {
  stopifnot(!anyNA(d$followup_state), all(d$followup_state %in%
    c("Conflicting evidence", "Death recorded", "Survival evidence", "Unresolved survival")))
  state <- ifelse(d$followup_state == "Conflicting evidence", "Conflicting evidence",
    ifelse(d$followup_state == "Death recorded", "Recorded death",
      ifelse(d$outcome_complete, "Both scores observed",
        ifelse(d$followup_state == "Survival evidence", "Alive, incomplete score pair", "Survival unresolved"))))
  categories <- c("Both scores observed", "Alive, incomplete score pair", "Recorded death",
                   "Survival unresolved", "Conflicting evidence")
  groups <- list(Training = d$training_eligible,
                 Complete_history = d$training_eligible & d$history_complete,
                 Incomplete_history = d$training_eligible & !d$history_complete)
  bind_rows(lapply(names(groups), function(nm) {
    keep <- groups[[nm]]
    counts <- as.integer(table(factor(state[keep], levels = categories)))
    data.frame(group = nm, category = categories, n = counts, denominator = sum(keep),
               percent = 100 * counts / sum(keep))
  }))
}

run_descriptive_outputs <- function(IN_COHORT = "output/01_cohort",
                                    IN_TRAJ = "output/02_trajectories",
                                    IN_MAIN = "output/03_main_models",
                                    IN_SENS = "output/04_sensitivity",
                                    OUT_DIR = "output/05_descriptive") {
  d <- readRDS(file.path(IN_COHORT, "charls_analytic_master.rds"))
  fc <- readRDS(file.path(IN_TRAJ, "full_cohort_with_trajectories_16717.rds"))
  train <- readRDS(file.path(IN_TRAJ, "landmark_training_with_trajectories_13851.rds"))
  lmk <- readRDS(file.path(IN_TRAJ, "landmark_primary_with_trajectories_7888.rds"))
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  csv <- function(x, name) write.csv(x, file.path(OUT_DIR, paste0(name, ".csv")), row.names = FALSE, na = "")
  groups <- list(Overall = rep(TRUE, nrow(fc)), `Longer-sleep` = fc$Sleep_Label == "Longer-sleep",
                 `Shorter-sleep` = fc$Sleep_Label == "Shorter-sleep")
  csv(observed_characteristics(fc, c(BLOCK_D, BLOCK_C), groups,
                               comparison = c("Shorter-sleep", "Longer-sleep")), "Table_1_observed")
  groups <- list(Included = train$ID %in% lmk$ID, Excluded = !train$ID %in% lmk$ID)
  csv(observed_characteristics(train, c("Sleep_Label", BLOCK_D, BLOCK_C,
                                        "c11", "c13", "c15", "c18", "c20", "n_sleep_exp"), groups),
      "Table_S3_included_excluded")
  # Both groups remain in the primary MI sample; this is not an exclusion table.
  stopifnot(!anyNA(lmk[setdiff(lmk_cols, c(BLOCK_D, BLOCK_C))]))
  cc <- complete.cases(lmk[c(BLOCK_D, BLOCK_C)])
  cc_data <- lmk
  cc_data$persistent_high_cesd <- factor(cc_data$persistent_high_cesd, levels = 0:1)
  csv(observed_characteristics(cc_data, c("Sleep_Label", BLOCK_D, BLOCK_C,
                                      "c11", "c13", "c15", "c18", "c20", "persistent_high_cesd"),
                               list(Complete = cc, Incomplete = !cc), c("Complete", "Incomplete")),
      "Table_S17_complete_incomplete_covariates")
  joint <- as.data.frame(table(Sleep = fc$Sleep_Label, Symptom = fc$Dep_Label), responseName = "n")
  joint$denominator <- ave(joint$n, joint$Sleep, FUN = sum)
  joint$percent_within_sleep <- 100 * joint$n / joint$denominator
  csv(joint, "Table_S4_class_cross_tabulation")
  csv(class_correspondence(fc), "class_correspondence_cluster_ci")
  csv(bind_rows(wave_counts(fc$n_sleep_5w, fc$Sleep_Label, "Sleep"),
                wave_counts(fc$n_cesd10_5w, fc$Dep_Label, "CES-D-10")), "Table_S4_observed_wave_counts")
  waves <- c(11, 13, 15, 18, 20)
  csv(bind_rows(lapply(c("s", "c"), function(prefix) {
    n <- colSums(!is.na(fc[paste0(prefix, waves)]))
    data.frame(domain = if (prefix == "s") "Sleep" else "CES-D-10", wave = 2000 + waves,
               n = unname(n), denominator = nrow(fc), percent = 100 * unname(n) / nrow(fc))
  })), "Table_S4_calendar_coverage")
  csv(bind_rows(posterior_quality(fc$Sleep_Label, fc$sleep_prob_assigned, "Full cohort sleep"),
                posterior_quality(fc$Dep_Label, fc$dep_prob_assigned, "Full cohort symptoms"),
                posterior_quality(train$Sleep_Label, train$sleep_prob_assigned, "Landmark training sleep")),
      "Table_S7_posterior_quality")
  csv(bind_rows(sleep_distribution(fc, paste0("s", waves), "Full cohort"),
                sleep_distribution(lmk, c("s11", "s13", "s15"), "Landmark regression")),
      "Table_S7_sleep_distribution")
  csv(bind_rows(missingness_summary(fc, "Full cohort"), missingness_summary(lmk, "Landmark primary")),
      "Table_S10_missingness")
  csv(cohort_flow(d), "Table_S13_eligibility")
  csv(followup_counts(d), "Table_S13_followup_status")
  primary <- readRDS(file.path(IN_MAIN, "primary_model_estimates.rds"))
  sensitivity <- readRDS(file.path(IN_SENS, "sensitivity_analysis_estimates.rds"))
  models <- list(`Full cohort Model 3` = primary$full_cohort$Model3,
                 `Primary landmark L5` = primary$landmark$L5)
  stems <- c("full_cohort_16717", "landmark_7888")
  mi_summaries <- lapply(seq_along(models), function(i) {
    mi <- readRDS(file.path(IN_MAIN, paste0("mids_", stems[i], ".rds")))
    completed <- readRDS(file.path(IN_MAIN, paste0("completed_", stems[i], ".rds")))
    mi_result_summary(models[[i]], mi, completed, names(models)[i])
  })
  csv(bind_rows(mi_summaries), "Table_S10B_mi_diagnostics")
  csv(bind_rows(lapply(names(models), function(nm) {
    cbind(analysis = nm, models[[nm]]$per_imputation)
  })), "Table_S10B_per_imputation")
  csv(estimate_rows(list(primary = primary, sensitivity = sensitivity)), "all_model_estimates")
  invisible(NULL)
}

if (sys.nframe() == 0L) run_descriptive_outputs()
