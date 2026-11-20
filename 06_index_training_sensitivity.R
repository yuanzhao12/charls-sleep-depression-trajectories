# Table S13D: restrict trajectory training to documented survival in 2015.
# Classification agreement only; no imputation or association-model refitting.
# Generated fits and outputs are private and must not be uploaded.

restrict_index_training_panel <- function(panel, training, primary_ids) {
  stopifnot(all(c("ID", "died2015") %in% names(training)),
            !anyNA(training$ID), !anyDuplicated(training$ID),
            setequal(panel$id_map$ID, training$ID))
  keep_ids <- training$ID[training$died2015 %in% 0]
  if (!all(primary_ids %in% keep_ids)) stop("The 2015-survival restriction excludes primary IDs.")
  # Subset the original training indices; do not rebuild or renumber the panel.
  panel$long <- panel$long[panel$long$ID %in% keep_ids, , drop = FALSE]
  panel$id_map <- panel$id_map[panel$id_map$ID %in% keep_ids, , drop = FALSE]
  stopifnot(nrow(panel$id_map) > 0L, setequal(panel$long$ID, panel$id_map$ID))
  panel
}

index_training_agreement <- function(training, primary, classes) {
  stopifnot(!anyNA(classes$ID), !anyDuplicated(classes$ID),
            all(classes$ID %in% training$ID), all(primary$ID %in% classes$ID))
  groups <- list(overlap_training = training[match(classes$ID, training$ID), , drop = FALSE],
                 primary_regression = primary)
  labels <- c("Longer-sleep", "Shorter-sleep")
  do.call(rbind, lapply(names(groups), function(population) {
    d <- groups[[population]]
    ref <- as.character(d$Sleep_Label)
    new <- as.character(classes$traj_label[match(d$ID, classes$ID)])
    stopifnot(nrow(d) > 0L, !anyDuplicated(d$ID),
              length(ref) == nrow(d), all(ref %in% labels), all(new %in% labels))
    tab <- table(factor(ref, levels = labels), factor(new, levels = labels))
    n <- sum(tab)
    pa <- sum(diag(tab)) / n
    pe <- sum(rowSums(tab) * colSums(tab)) / n^2
    data.frame(population = population, n = n, agreement = pa,
               kappa = if (pe < 1) (pa - pe) / (1 - pe) else NA_real_,
               changed = sum(new != ref), primary_short = sum(ref == "Shorter-sleep"),
               index2015_short = sum(new == "Shorter-sleep"))
  }))
}

run_index_training_sensitivity <- function(in_dir = "output/02_trajectories",
                                           out_dir = "output/06_index_training_sensitivity",
                                           trajectory_script = "02_trajectory_models.R") {
  trajectory <- new.env(parent = environment())
  sys.source(trajectory_script, envir = trajectory)
  training <- readRDS(file.path(in_dir, "landmark_training_with_trajectories_13851.rds"))
  primary <- readRDS(file.path(in_dir, "landmark_primary_with_trajectories_7888.rds"))
  stopifnot(nrow(training) == 13851L, nrow(primary) == 7888L,
            all(c("ID", "Sleep_Label") %in% names(primary)),
            all(c("ID", "Sleep_Label", "died2015", "s11", "s13", "s15") %in% names(training)),
            !anyNA(primary$ID), !anyDuplicated(primary$ID), all(primary$ID %in% training$ID),
            all(training$Sleep_Label %in% c("Longer-sleep", "Shorter-sleep")),
            all(as.character(primary$Sleep_Label) ==
                as.character(training$Sleep_Label[match(primary$ID, training$ID)])))

  panel <- trajectory$build_long_panel(training, c("s11", "s13", "s15"),
                                        c(2011L, 2013L, 2015L), min_obs = 2L)
  panel <- restrict_index_training_panel(panel, training, primary$ID)
  stopifnot(nrow(panel$id_map) == 13345L)
  result <- trajectory$fit_trajectory_grid(panel, file.path(out_dir, "fit"),
                                           max_k = 2L, seeds = 2025:2027, workers = 24L)
  classes <- trajectory$assign_neutral_labels(result$selected, panel,
    lower_label = "Shorter-sleep", higher_label = "Longer-sleep", ref_label = "Longer-sleep")
  agreement <- index_training_agreement(training, primary, classes)
  utils::write.csv(agreement, file.path(out_dir, "TRAINING_AGREEMENT.csv"), row.names = FALSE, na = "")
  invisible(list(agreement = agreement, classes = classes, grid_table = result$grid_table))
}

if (sys.nframe() == 0L) run_index_training_sensitivity()
