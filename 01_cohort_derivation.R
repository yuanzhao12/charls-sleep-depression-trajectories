# Cohort linkage, baseline covariates, observed sleep and symptom scores.

library(haven)
library(dplyr)

# Canonical 12-digit CHARLS participant identifier
norm_charls_id <- function(x, allow_missing = FALSE) {
  x <- trimws(as.character(x))
  idx <- !is.na(x) & nchar(x) == 11L
  x[idx] <- paste0(substr(x[idx], 1L, 9L), "0", substr(x[idx], 10L, 11L))
  bad <- is.na(x) | !grepl("^[0-9]{12}$", x)
  if (!allow_missing && any(bad)) stop("Missing or invalid participant identifier")
  x[bad] <- NA_character_
  stopifnot(!anyDuplicated(x[!is.na(x)]))
  x
}

read_charls_dta <- function(root, rel_path, cols = NULL, allow_missing_id = FALSE) {
  full_path <- file.path(root, rel_path)
  df <- if (is.null(cols)) {
    read_dta(full_path)
  } else {
    read_dta(full_path, col_select = all_of(unique(c("ID", cols))))
  }
  if ("ID" %in% names(df)) df$ID <- norm_charls_id(df$ID, allow_missing_id)
  df
}

# Nighttime sleep duration restricted to plausible range (2 to 12 hours/night)
clean_sleep_hours <- function(x) {
  v <- as.numeric(x)
  v[!is.finite(v) | v < 2 | v > 12] <- NA_real_
  v
}

# CES-D-10 (0-30) and CES-D-9 excluding restless sleep (0-27) with prorated scoring
score_cesd <- function(df, item_cols, drop_sleep = FALSE, min_items = 8L) {
  stopifnot(length(item_cols) == 10L, min_items >= 1L, min_items <= if (drop_sleep) 9L else 10L)
  mat <- vapply(df[item_cols], as.numeric, numeric(nrow(df)))
  mat[!mat %in% 1:4] <- NA_real_
  mat <- mat - 1

  # Reverse-score positive affect items: item 5 (hopeful) and item 8 (happy)
  mat[, c(5L, 8L)] <- 3 - mat[, c(5L, 8L)]

  if (drop_sleep) {
    mat <- mat[, -7L, drop = FALSE]
    n_valid <- rowSums(!is.na(mat))
    score <- ifelse(n_valid >= min_items, rowSums(mat, na.rm = TRUE) * (9 / n_valid), NA_real_)
    return(list(score = score, n_valid = n_valid))
  }

  n_valid <- rowSums(!is.na(mat))
  score <- ifelse(n_valid >= min_items, rowSums(mat, na.rm = TRUE) * (10 / n_valid), NA_real_)
  list(score = score, n_valid = n_valid)
}

# Blood pressure readings 2 and 3 with measurement-failure (993/999) and pulse-pressure rules
clean_bp_means <- function(s2, s3, d2, d3) {
  s2 <- as.numeric(s2); s3 <- as.numeric(s3)
  d2 <- as.numeric(d2); d3 <- as.numeric(d3)
  failed2 <- s2 %in% c(993, 999)
  failed3 <- s3 %in% c(993, 999)

  s2[!is.finite(s2) | s2 < 60 | s2 > 300] <- NA_real_
  s3[!is.finite(s3) | s3 < 60 | s3 > 300] <- NA_real_
  d2[failed2 | !is.finite(d2) | d2 < 30 | d2 > 150] <- NA_real_
  d3[failed3 | !is.finite(d3) | d3 < 30 | d3 > 150] <- NA_real_

  pp2 <- !is.na(s2) & !is.na(d2) & s2 <= d2
  pp3 <- !is.na(s3) & !is.na(d3) & s3 <= d3
  s2[pp2] <- d2[pp2] <- NA_real_
  s3[pp3] <- d3[pp3] <- NA_real_

  list(SBP = (s2 + s3) / 2, DBP = (d2 + d3) / 2)
}

yn_logical <- function(x) c(TRUE, FALSE)[match(as.numeric(x), 1:2)]
yn_factor  <- function(x) factor(x, levels = c(FALSE, TRUE), labels = c("No", "Yes"))

run_cohort_derivation <- function() {
  RAW_ROOT <- Sys.getenv("CHARLS_RAW_ROOT", unset = "data/raw_charls")
  STATUS_ROOT <- Sys.getenv("CHARLS_STATUS_ROOT", unset = RAW_ROOT)
  OUT_DIR <- "output/01_cohort"
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
  bg11     <- read_charls_dta(RAW_ROOT, "2011/demographic_background.dta",
                              c("householdID", "communityID", "ba002_1", "ba004", "rgender", "bd001", "be001"))
  hsf11    <- read_charls_dta(RAW_ROOT, "2011/health_status_and_functioning.dta",
                              c("householdID", "communityID", "da049", sprintf("dc0%02d", 9:18),
                                "da059", "da061", "da067", "da007_1_", "da008_1_", "da007_2_",
                                "da007_3_", "da007_7_", "da007_8_", "da041", sprintf("db0%02d", 10:15)))
  bio11    <- read_charls_dta(RAW_ROOT, "2011/biomarkers.dta",
                              c("householdID", "communityID", "qa007", "qa011", "qa008", "qa012", "ql002", "qi002"))
  blood11  <- read_charls_dta(RAW_ROOT, "2011/Blood_20140429.dta",
                              c("newtg", "newhdl", "newcho", "newldl", "newglu", "newhba1c", "qc1_va003"))
  roster11 <- read_charls_dta(RAW_ROOT, "2011/household_roster.dta",
                              c("householdID", "communityID", "a001"))
  weight11 <- read_charls_dta(RAW_ROOT, "2011/weight.dta", c("communityID", "ind_weight_ad2"),
                              allow_missing_id = TRUE)
  stopifnot(!any(is.na(weight11$ID) & is.finite(weight11$ind_weight_ad2) & weight11$ind_weight_ad2 > 0))
  survey_design <- as.data.frame(weight11[!is.na(weight11$ID) &
    is.finite(weight11$ind_weight_ad2) & weight11$ind_weight_ad2 > 0, ])
  survey_design$communityID <- as.character(survey_design$communityID)
  survey_design$base_weight <- as.numeric(survey_design$ind_weight_ad2)
  survey_design <- survey_design[c("ID", "communityID", "base_weight")]
  stopifnot(!anyNA(survey_design), !anyDuplicated(survey_design$ID),
            length(unique(survey_design$communityID)) == 450L)
  
  baseline_ids <- sort(Reduce(union, list(bg11$ID, hsf11$ID, bio11$ID, blood11$ID, roster11$ID)))
  d <- data.frame(ID = baseline_ids, stringsAsFactors = FALSE)
  
  gi <- match(d$ID, bg11$ID)
  hi <- match(d$ID, hsf11$ID)
  mi <- match(d$ID, bio11$ID)
  bi <- match(d$ID, blood11$ID)
  wi <- match(d$ID, weight11$ID)
  stopifnot(!anyDuplicated(roster11$householdID))
  hh <- match(bg11$householdID[gi], roster11$householdID)
  
  d$householdID    <- as.character(bg11$householdID[gi])
  d$communityID    <- as.character(bg11$communityID[gi])
  d$ind_weight_ad2 <- as.numeric(weight11$ind_weight_ad2[wi])
  
  # Block D: Sociodemographic and lifestyle factors
  d$age <- ifelse(is.na(bg11$ba002_1[gi]), as.numeric(bg11$ba004[gi]), 2011 - as.numeric(bg11$ba002_1[gi]))
  sex_v <- as.numeric(bg11$rgender[gi])
  d$gender <- factor(ifelse(sex_v %in% 1:2, ifelse(sex_v == 1, "Male", "Female"), NA_character_),
                     levels = c("Male", "Female"))
  
  edu_levels <- c("Below primary", "Primary", "Middle school or above")
  edu_groups <- c(rep(1L, 2L), rep(2L, 2L), rep(3L, 7L))
  d$education <- factor(edu_levels[edu_groups[match(as.numeric(bg11$bd001[gi]), 1:11)]], levels = edu_levels)
  
  mar_v <- as.numeric(bg11$be001[gi])
  d$marital <- factor(c(rep("Married", 2L), rep("Unmarried", 4L))[match(mar_v, 1:6)],
                      levels = c("Married", "Unmarried"))
  
  area_v <- as.numeric(roster11$a001[hh])
  d$residence <- factor(c("Rural", "Urban")[match(area_v, 1:2)], levels = c("Rural", "Urban"))
  
  smk_ever <- as.numeric(hsf11$da059[hi])
  smk_now  <- as.numeric(hsf11$da061[hi])
  d$smoking <- factor(
    ifelse(smk_ever %in% 2, "Never smoker",
      ifelse(smk_ever %in% 1 & smk_now %in% 2, "Former smoker",
        ifelse(smk_ever %in% 1 & smk_now %in% 1, "Current smoker", NA_character_))),
    levels = c("Never smoker", "Former smoker", "Current smoker")
  )
  
  drink_levels <- c("None", "Less than once a month", "More than once a month")
  d$drinking <- factor(rev(drink_levels)[match(as.numeric(hsf11$da067[hi]), 1:3)], levels = drink_levels)
  
  # Block C: BMI and cardiometabolic conditions (TRUE overrides unknown; FALSE requires observed negative)
  kg <- as.numeric(bio11$ql002[mi])
  cm <- as.numeric(bio11$qi002[mi])
  d$BMI <- ifelse(kg > 30 & kg < 150, kg, NA_real_) / (ifelse(cm > 100 & cm < 200, cm, NA_real_) / 100)^2
  
  bp <- clean_bp_means(bio11$qa007[mi], bio11$qa011[mi], bio11$qa008[mi], bio11$qa012[mi])
  d$SBP <- bp$SBP
  d$DBP <- bp$DBP
  
  clean_pos <- function(v) { z <- as.numeric(v); replace(z, !is.finite(z) | z <= 0, NA_real_) }
  d$TG      <- clean_pos(blood11$newtg[bi])
  d$HDL     <- clean_pos(blood11$newhdl[bi])
  d$TC      <- clean_pos(blood11$newcho[bi])
  d$LDL     <- clean_pos(blood11$newldl[bi])
  d$Glucose <- clean_pos(blood11$newglu[bi])
  d$HbA1c   <- clean_pos(blood11$newhba1c[bi])
  fast_v    <- as.numeric(blood11$qc1_va003[bi])
  d$FPG     <- replace(d$Glucose, !fast_v %in% 1, NA_real_)
  
  htn_sr <- yn_logical(hsf11$da007_1_[hi]) | yn_logical(hsf11$da008_1_[hi])
  dm_sr  <- yn_logical(hsf11$da007_3_[hi])
  lip_sr <- yn_logical(hsf11$da007_2_[hi])
  
  d$HTN          <- yn_factor(htn_sr | d$SBP >= 140 | d$DBP >= 90)
  d$diabetes     <- yn_factor(dm_sr  | d$HbA1c >= 6.5 | d$FPG >= 126)
  d$dyslipidemia <- yn_factor(lip_sr | d$TC >= 240 | d$LDL >= 160 | d$HDL < 40 | d$TG >= 200)
  d$heart_hist   <- yn_factor(yn_logical(hsf11$da007_7_[hi]))
  d$stroke_hist  <- yn_factor(yn_logical(hsf11$da007_8_[hi]))
  
  # Broader health sensitivity covariates: baseline pain and 6-item basic ADL difficulty
  d$pain <- yn_factor(yn_logical(hsf11$da041[hi]))
  adl_cols <- sprintf("db0%02d", 10:15)
  adl_mat  <- vapply(hsf11[hi, adl_cols], function(v) {
    z <- as.numeric(v); replace(z, !z %in% 1:4, NA_real_)
  }, numeric(nrow(d)))
  adl_obs6 <- rowSums(!is.na(adl_mat)) == 6L
  d$adl_limit <- yn_factor(ifelse(adl_obs6, rowSums(adl_mat > 1) > 0L, NA))
  
  # Five-wave sleep duration and CES-D symptom scoring (2011, 2013, 2015, 2018, 2020)
  wave_specs <- list(
    `2011` = list(hsf = "2011/health_status_and_functioning.dta", cog = NULL,                 sleep = "da049", items = sprintf("dc0%02d", 9:18)),
    `2013` = list(hsf = "2013/Health_Status_and_Functioning.dta", cog = NULL,                 sleep = "da049", items = sprintf("dc0%02d", 9:18)),
    `2015` = list(hsf = "2015/Health_Status_and_Functioning.dta", cog = NULL,                 sleep = "da049", items = sprintf("dc0%02d", 9:18)),
    `2018` = list(hsf = "2018/Health_Status_and_Functioning.dta", cog = "2018/Cognition.dta", sleep = "da049", items = sprintf("dc0%02d", 9:18)),
    `2020` = list(hsf = "2020/Health_Status_and_Functioning.dta", cog = NULL,                 sleep = "da030", items = sprintf("dc0%02d", 16:25))
  )
  
  for (yr in names(wave_specs)) {
    sp  <- wave_specs[[yr]]
    suf <- substr(yr, 3L, 4L)
    s_df <- read_charls_dta(RAW_ROOT, sp$hsf, sp$sleep)
    c_df <- if (is.null(sp$cog)) read_charls_dta(RAW_ROOT, sp$hsf, sp$items) else read_charls_dta(RAW_ROOT, sp$cog, sp$items)
  
    sc10 <- score_cesd(c_df, sp$items, drop_sleep = FALSE, min_items = 8L)
    sc9  <- score_cesd(c_df, sp$items, drop_sleep = TRUE,  min_items = 8L)
    sc9s <- score_cesd(c_df, sp$items, drop_sleep = TRUE,  min_items = 9L)
  
    ix_c <- match(d$ID, c_df$ID)
    ix_s <- match(d$ID, s_df$ID)
  
    d[[paste0("s", suf)]]        <- clean_sleep_hours(s_df[[sp$sleep]][ix_s])
    d[[paste0("c", suf)]]        <- sc10$score[ix_c]
    d[[paste0("c9_", suf)]]      <- sc9$score[ix_c]
    d[[paste0("c9strict_", suf)]] <- sc9s$score[ix_c]
  }
  
  # Vital status reconciliation at 2015 index and 2018-2020 follow-up waves
  ex13 <- read_charls_dta(RAW_ROOT, "2013/Exit_Interview.dta", character())
  va13 <- read_charls_dta(RAW_ROOT, "2013/Verbal_Autopsy.dta", character())
  si15 <- read_charls_dta(STATUS_ROOT, "2015/Sample_Infor.dta", c("died", "iyear", "imonth"))
  si18 <- read_charls_dta(STATUS_ROOT, "2018/Sample_Infor.dta", c("died", "iyear", "imonth"))
  si20 <- read_charls_dta(STATUS_ROOT, "2020/Sample_Infor.dta", c("died", "iyear", "imonth"))
  ex20 <- read_charls_dta(STATUS_ROOT, "2020/Exit_Module.dta",  c("exb001_1", "exb001_2", "exb001_3"))
  
  d$died2015 <- as.numeric(si15$died[match(d$ID, si15$ID)])
  d$died2018 <- as.numeric(si18$died[match(d$ID, si18$ID)])
  d$died2020 <- as.numeric(si20$died[match(d$ID, si20$ID)])
  
  i_exit     <- match(d$ID, ex20$ID)
  dy         <- as.integer(ex20$exb001_1[i_exit])
  dm         <- as.integer(ex20$exb001_2[i_exit])
  early_date <- is.finite(dy) & dy >= 1900L & dy < 2015L
  early_rec  <- d$ID %in% union(ex13$ID, va13$ID)
  
  w3_alive    <- d$died2015 %in% 0L
  later_alive <- d$died2018 %in% 0L | d$died2020 %in% 0L
  w3_score    <- !is.na(d$c15)
  
  conflict_idx <- ((d$died2015 %in% 1L | early_date | early_rec) & (w3_alive | later_alive | w3_score)) |
                  (d$died2018 %in% 1L & d$died2020 %in% 0L)
  dead_by_idx  <- d$died2015 %in% 1L | early_date | early_rec
  alive_idx    <- w3_alive | later_alive | w3_score
  
  d$vital_index <- ifelse(conflict_idx, "CONFLICT",
                     ifelse(dead_by_idx, "DEAD_BY_INDEX",
                       ifelse(alive_idx, "ALIVE_EVIDENCE", "UNKNOWN")))
  
  # Post-2015 vital classification for 2018 & 2020 outcome window
  earlier_death <- function(dy, dm, iy, im) (dy < iy) %in% TRUE | (dy == iy & dm < im) %in% TRUE
  status_waves  <- list(`2015` = si15, `2018` = si18, `2020` = si20)
  date_conflicts <- matrix(FALSE, nrow(d), 3L, dimnames = list(NULL, names(status_waves)))
  
  for (wy in names(status_waves)) {
    sw <- status_waves[[wy]]
    ix <- match(d$ID, sw$ID)
    pos <- (sw$died[ix] %in% 0L) | !is.na(d[[paste0("c", substr(wy, 3L, 4L))]])
    iy  <- as.integer(sw$iyear[ix])
    im  <- as.integer(sw$imonth[ix])
    iy[pos & is.na(iy)] <- as.integer(wy)
    date_conflicts[, wy] <- pos & earlier_death(dy, dm, iy, im)
  }
  
  obs18 <- !is.na(d$c18)
  obs20 <- !is.na(d$c20)
  ev18  <- d$died2018 %in% 0L | obs18
  ev20  <- d$died2020 %in% 0L | obs20
  dead18 <- d$died2018 %in% 1L
  dead20 <- d$died2020 %in% 1L | !is.na(i_exit)
  conflict_post <- (dead18 & (ev18 | ev20)) | (dead20 & ev20) | (rowSums(date_conflicts) > 0L)
  dead_post     <- dead18 | dead20
  d$alive_both_2020 <- !conflict_post & !dead_post & ev20
  d$followup_state <- ifelse(conflict_post, "Conflicting evidence",
    ifelse(dead_post, "Death recorded", ifelse(ev20, "Survival evidence", "Unresolved survival")))
  
  # Analytical cohort eligibility flags and outcomes
  sleep_5w  <- c("s11", "s13", "s15", "s18", "s20")
  cesd10_5w <- c("c11", "c13", "c15", "c18", "c20")
  cesd9_5w  <- c("c9_11", "c9_13", "c9_15", "c9_18", "c9_20")
  cesd9s_5w <- c("c9strict_11", "c9strict_13", "c9strict_15", "c9strict_18", "c9strict_20")
  
  d$n_sleep_5w       <- rowSums(!is.na(d[sleep_5w]))
  d$n_cesd10_5w      <- rowSums(!is.na(d[cesd10_5w]))
  d$n_cesd9_5w       <- rowSums(!is.na(d[cesd9_5w]))
  d$n_sleep_exp      <- rowSums(!is.na(d[c("s11", "s13", "s15")]))
  d$mean_sleep_11_15 <- rowMeans(d[c("s11", "s13", "s15")], na.rm = TRUE)
  
  # Full-cohort eligibility (n = 16,717; paired 9-item n = 16,693)
  d$eligible_45plus      <- !is.na(d$age) & d$age >= 45 & !is.na(d$gender)
  d$in_full_cohort       <- d$eligible_45plus & d$n_sleep_5w >= 1L & d$n_cesd10_5w >= 1L
  d$in_full_paired9      <- d$in_full_cohort  & d$n_cesd9_5w >= 1L
  
  # Landmark eligibility (sleep-eligible n = 14,426; training n = 13,851; primary n = 7,888)
  d$age_sleep_eligible   <- d$eligible_45plus & d$n_sleep_exp >= 2L
  d$training_eligible    <- d$age_sleep_eligible & d$vital_index == "ALIVE_EVIDENCE"
  d$history_complete     <- complete.cases(d[c("c11", "c13", "c15")])
  d$outcome_complete     <- obs18 & obs20
  d$common_eligible      <- d$training_eligible & d$history_complete & d$outcome_complete
  d$strict_eligible      <- d$common_eligible & (d$c11 < 10 & d$c13 < 10 & d$c15 < 10)
  d$paired8_eligible     <- d$common_eligible & complete.cases(d[cesd9_5w])
  d$paired9_eligible     <- d$common_eligible & complete.cases(d[cesd9s_5w])
  d$survivor_resp_target <- d$training_eligible & d$history_complete & d$alive_both_2020
  
  # Landmark outcomes: binary persistent elevation (2018 & 2020 CES-D-10 >= 10) and continuous two-visit means
  d$persistent_high_cesd <- ifelse(d$outcome_complete, as.integer(d$c18 >= 10 & d$c20 >= 10), NA_integer_)
  d$mean10_1820          <- ifelse(d$outcome_complete, (d$c18 + d$c20) / 2, NA_real_)
  d$mean9_1820           <- ifelse(complete.cases(d[c("c9_18", "c9_20")]), (d$c9_18 + d$c9_20) / 2, NA_real_)
  d$mean9_strict_1820    <- ifelse(complete.cases(d[c("c9strict_18", "c9strict_20")]),
                                   (d$c9strict_18 + d$c9strict_20) / 2, NA_real_)
  
  # Published counts are checks on the data release, not selection targets.
  counts <- data.frame(
    sample = c("Baseline records", "Full cohort", "Full-cohort paired nine-item", "Landmark sleep-eligible",
               "Landmark training", "Landmark primary", "Below-threshold subset", "Paired eight-of-nine",
               "Complete nine-item", "Survivor response target"),
    n = c(nrow(d), sum(d$in_full_cohort), sum(d$in_full_paired9), sum(d$age_sleep_eligible),
          sum(d$training_eligible), sum(d$common_eligible), sum(d$strict_eligible),
          sum(d$paired8_eligible), sum(d$paired9_eligible), sum(d$survivor_resp_target)))
  expected <- c(17713L, 16717L, 16693L, 14426L, 13851L, 7888L, 3457L, 7446L, 6224L, 9297L)
  write.csv(counts, file.path(OUT_DIR, "cohort_counts.csv"), row.names = FALSE)
  if (!identical(as.integer(counts$n), expected)) stop("Cohort counts differ from the analysed release; inspect cohort_counts.csv.")
  stopifnot(sum(d$persistent_high_cesd[d$common_eligible]) == 2051L,
            sum(d$persistent_high_cesd[d$paired8_eligible]) == 1901L,
            sum(d$persistent_high_cesd[d$paired9_eligible]) == 1513L,
            !anyNA(d$communityID[d$in_full_cohort]))
  
  saveRDS(d, file.path(OUT_DIR, "charls_analytic_master.rds"))
  saveRDS(d[d$in_full_cohort, ],    file.path(OUT_DIR, "cohort_full_16717.rds"))
  saveRDS(d[d$training_eligible, ], file.path(OUT_DIR, "cohort_landmark_training_13851.rds"))
  saveRDS(d[d$common_eligible, ],   file.path(OUT_DIR, "cohort_landmark_primary_7888.rds"))
  saveRDS(survey_design, file.path(OUT_DIR, "baseline_survey_design.rds"))
  
  invisible(d)
}

if (sys.nframe() == 0L) run_cohort_derivation()
