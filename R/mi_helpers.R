# Baseline-covariate imputation. Longitudinal scores and class assignments stay observed.

check_fcs_design <- function(y, ry, x, wy, equation) {
  stopifnot(is.matrix(x), nrow(x) == equation$n,
            identical(colnames(x), equation$columns), all(is.finite(x)),
            identical(which(unname(wy)), equation$missing_rows),
            identical(which(!unname(ry)), equation$missing_rows))
  X <- cbind(1, x[ry, , drop = FALSE])
  stopifnot(nrow(X) > ncol(X), qr(X)$rank == ncol(X), length(unique(y[ry])) > 1L)
  if (is.factor(y)) stopifnot(all(tabulate(as.integer(y[ry]), nlevels(y)) > 0L))
}

# Inspect convergence without changing mice's multinomial fit or random draws.
checked_categorical_draw <- function(method, y, ry, x, wy, ...) {
  fun <- switch(method, polyreg = mice::mice.impute.polyreg, polr = mice::mice.impute.polr)
  expr <- body(fun)
  at <- which(vapply(as.list(expr)[-1L], function(z) {
    is.call(z) && length(z) == 3L && identical(z[[1L]], as.name("<-")) &&
      identical(z[[2L]], quote(fit))
  }, logical(1))) + 1L
  stopifnot(length(at) == 1L)
  call <- expr[[at]][[3L]]
  backend <- if (method == "polr") quote(MASS::polr) else quote(nnet::multinom)
  if (method == "polr") {
    stopifnot(identical(call[[2L]][[2L]][[1L]], backend))
    call[[2L]][[2L]][[1L]] <- quote(.checked_optimizer)
  } else {
    stopifnot(identical(call[[1L]], backend))
    call[[1L]] <- quote(.checked_optimizer)
  }
  status <- NULL
  scope <- new.env(parent = environment(fun))
  scope$.checked_optimizer <- function(...) {
    original <- sys.call()
    original[[1L]] <- backend
    fit <- eval(original, envir = parent.frame())
    status <<- fit$convergence
    fit
  }
  expr[[at]][[3L]] <- call
  body(fun) <- expr
  environment(fun) <- scope
  draws <- fun(y = y, ry = ry, x = x, wy = wy, ...)
  if (!isTRUE(status == 0)) stop("Categorical imputation did not converge: ", method)
  draws
}

mice.impute.checked_pmm <- function(y, ry, x, wy, equation, ...) {
  check_fcs_design(y, ry, x, wy, equation)
  mice::mice.impute.pmm(y, ry, x, wy, rank.adjust = FALSE)
}
mice.impute.checked_logreg <- function(y, ry, x, wy, equation, ...) {
  check_fcs_design(y, ry, x, wy, equation)
  mice::mice.impute.logreg(y, ry, x, wy)
}
mice.impute.checked_polyreg <- function(y, ry, x, wy, equation, ...) {
  check_fcs_design(y, ry, x, wy, equation)
  checked_categorical_draw("polyreg", y, ry, x, wy)
}
mice.impute.checked_polr <- function(y, ry, x, wy, equation, ...) {
  check_fcs_design(y, ry, x, wy, equation)
  checked_categorical_draw("polr", y, ry, x, wy, polr.to.loggedEvents = TRUE)
}

run_covariate_mice <- function(dat, active_predictors, m = 100L, maxit = 20L,
                               seed = 2025L, profile = c("standard", "primary", "sensitivity")) {
  profile <- match.arg(profile)
  stopifnot(as.character(utils::packageVersion("mice")) == "3.19.0")
  vars <- names(dat)
  stopifnot(!anyDuplicated(vars), all(active_predictors %in% vars))
  targets <- intersect(IMPUTABLE, vars)
  targets <- targets[vapply(dat[targets], anyNA, logical(1))]
  stopifnot(length(targets) > 0L, !anyNA(dat[setdiff(vars, targets)]))
  P <- matrix(0L, length(vars), length(vars), dimnames = list(vars, vars))
  for (v in targets) P[v, setdiff(active_predictors, v)] <- 1L
  method <- setNames(rep("", length(vars)), vars)
  for (v in targets) {
    x <- dat[[v]]
    stopifnot(!all(is.na(x)), is.numeric(x) || is.factor(x))
    method[v] <- if (is.numeric(x)) "pmm" else if (nlevels(x) == 2L) "logreg" else
      if (is.ordered(x)) "polr" else "polyreg"
  }
  engine_method <- method
  args <- list(data = dat, m = m, maxit = maxit, seed = seed,
               predictorMatrix = P, visitSequence = targets, where = is.na(dat), printFlag = FALSE)
  if (profile == "primary") {
    engine_method[targets] <- paste0("checked_", method[targets])
    args$blots <- setNames(lapply(vars, function(v) {
      if (!v %in% targets) return(list())
      form <- reformulate(names(which(P[v, ] != 0L)))
      mf <- model.frame(form, data = dat, na.action = na.pass)
      cols <- colnames(model.matrix(form, mf))[-1L]
      list(equation = list(n = nrow(dat), columns = cols, missing_rows = which(is.na(dat[[v]]))))
    }), vars)
    args$remove.constant <- FALSE
    args$remove.collinear <- FALSE
    args$eps <- 0
  } else if (profile == "standard") {
    args$nnet.maxit <- 200L
  } else {
    args$polr.to.loggedEvents <- TRUE
  }
  args$method <- engine_method
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  mi <- do.call(mice::mice, args)
  stopifnot(identical(mi$method, engine_method), all(mi$predictorMatrix == P),
            identical(unname(mi$visitSequence), targets), NROW(mi$loggedEvents) == 0L,
            mi$m == m, mi$iteration == maxit, all(mi$where == is.na(dat)))
  for (j in seq_len(m)) {
    completed <- mice::complete(mi, j)
    stopifnot(identical(names(completed), vars), nrow(completed) == nrow(dat), !anyNA(completed))
    for (v in vars) {
      observed <- !is.na(dat[[v]])
      stopifnot(identical(completed[[v]][observed], dat[[v]][observed]))
    }
  }
  mi
}

write_mi_diagnostics <- function(mi, out_dir, prefix) {
  methods <- data.frame(variable = names(mi$method), method = unname(mi$method),
                        missing = unname(mi$nmis), visit_order = match(names(mi$method), mi$visitSequence))
  write.csv(methods, file.path(out_dir, paste0(prefix, "_methods.csv")), row.names = FALSE)
  write.csv(data.frame(target = rownames(mi$predictorMatrix), mi$predictorMatrix, check.names = FALSE),
            file.path(out_dir, paste0(prefix, "_predictors.csv")), row.names = FALSE)
  chain <- as.data.frame.table(mi$chainMean, responseName = "chain_mean", stringsAsFactors = FALSE)
  names(chain)[1:3] <- c("variable", "iteration", "imputation")
  chain$chain_variance <- as.vector(mi$chainVar)
  chain <- chain[chain$variable %in% mi$visitSequence, ]
  write.csv(chain, file.path(out_dir, paste0(prefix, "_chains.csv")), row.names = FALSE)
  events <- mi$loggedEvents
  if (is.null(events)) events <- data.frame(iteration = integer(), imputation = integer(), detail = character())
  write.csv(events, file.path(out_dir, paste0(prefix, "_events.csv")), row.names = FALSE)
}
