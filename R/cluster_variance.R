# Named covariance validation and MI pooling; no fitting or I/O on import.
charls_check_vcov <- function(v, terms) {
  if (!is.matrix(v) || !identical(dim(v), rep(length(terms), 2L)) ||
      !identical(rownames(v), terms) || !identical(colnames(v), terms) ||
      any(!is.finite(v)) || max(abs(v - t(v))) > 1e-8) {
    stop("invalid covariance dimensions, names, order or values", call. = FALSE)
  }
  if (min(eigen(v, symmetric = TRUE, only.values = TRUE)$values) < -1e-8) {
    stop("covariance matrix is not positive semidefinite", call. = FALSE)
  }
  v
}

charls_cluster_vcov <- function(fit, cluster) {
  ids <- rownames(stats::model.frame(fit))
  if (is.null(names(cluster)) || anyDuplicated(names(cluster)) ||
      !setequal(names(cluster), ids) || anyNA(cluster) || any(!nzchar(as.character(cluster)))) {
    stop("cluster must be a complete named vector keyed by fitted subject IDs", call. = FALSE)
  }
  cl <- as.character(cluster[ids])
  if (length(unique(cl)) < 2L) stop("at least two clusters required", call. = FALSE)
  v <- sandwich::vcovCL(fit, cluster = factor(cl), type = "HC1", cadjust = TRUE)
  v <- charls_check_vcov(v, names(stats::coef(fit)))
  list(vcov = v, n_clusters = length(unique(cl)), cluster_sizes = table(cl),
       fitted_ids = ids, type = "HC1", cadjust = TRUE)
}

charls_rubin_pool <- function(est, var, dfcom = Inf,
                              scale = c("log_odds", "identity", "log_ratio"),
                              conf.level = 0.95) {
  scale <- match.arg(scale)
  m <- length(est)
  if (!is.list(est) || !is.list(var) || m < 2L || length(var) != m) {
    stop("pooling requires matching lists of at least two imputations", call. = FALSE)
  }
  if (length(dfcom) != 1L || is.na(dfcom) || dfcom <= 0 ||
      length(conf.level) != 1L || !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1) {
    stop("invalid complete-data degrees of freedom or confidence level", call. = FALSE)
  }
  terms <- names(est[[1]])
  if (is.null(terms) || !length(terms) || anyDuplicated(terms) || any(!nzchar(terms))) {
    stop("coefficient names must be nonempty and unique", call. = FALSE)
  }
  p <- length(terms)
  covariance <- lapply(seq_len(m), function(i) {
    if (!is.numeric(est[[i]]) || !identical(names(est[[i]]), terms) || any(!is.finite(est[[i]]))) {
      stop("different coefficient names/order or nonfinite estimates", call. = FALSE)
    }
    v <- var[[i]]
    if (!is.matrix(v)) {
      if (!is.numeric(v) || !identical(names(v), terms) || any(!is.finite(v)) || any(v < 0)) {
        stop("variance vector names/order or values invalid", call. = FALSE)
      }
      v <- diag(as.numeric(v), nrow = p, ncol = p)
      dimnames(v) <- list(terms, terms)
    }
    charls_check_vcov(v, terms)
  })
  E <- do.call(rbind, est)
  Ubar <- Reduce(`+`, covariance) / m
  B <- matrix(stats::cov(E), p, p, dimnames = list(terms, terms))
  Tmat <- Ubar + (1 + 1 / m) * B
  rows <- lapply(seq_len(p), function(j) {
    q <- E[, j]
    u <- vapply(covariance, function(v) v[j, j], numeric(1))
    if (all(u == 0) && stats::var(q) == 0) {
      qbar <- mean(q); se <- 0; df <- dfcom; fmi <- 0; lambda <- 0
      lower <- upper <- qbar; statistic <- pvalue <- NA_real_
    } else {
      # Reuse the tested Barnard-Rubin implementation instead of re-deriving df.
      ref <- mice::pool.scalar(q, u, n = dfcom + 1, k = 1, rule = "rubin1987")
      qbar <- ref$qbar; se <- sqrt(ref$t); df <- ref$df; fmi <- ref$fmi
      lambda <- (1 + 1 / m) * ref$b / ref$t
      critical <- stats::qt((1 + conf.level) / 2, df)
      lower <- qbar - critical * se; upper <- qbar + critical * se
      statistic <- qbar / se
      pvalue <- 2 * stats::pt(-abs(statistic), df = df)
    }
    data.frame(term = terms[j], estimate = qbar, se = se, lower = lower,
               upper = upper, statistic = statistic, p = pvalue, df = df,
               dfcom = dfcom, lambda = lambda, fmi = fmi, scale = scale,
               zero_variance = se == 0, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (scale == "log_odds") {
    out$OR <- exp(out$estimate); out$OR_lower <- exp(out$lower); out$OR_upper <- exp(out$upper)
  } else if (scale == "log_ratio") {
    out$RR <- exp(out$estimate); out$RR_lower <- exp(out$lower); out$RR_upper <- exp(out$upper)
  }
  list(table = out, m = m, Ubar = Ubar, B = B, T = Tmat,
       lambda = out$lambda, scale = scale, confidence_level = conf.level,
       note = "Named within-imputation covariance retained; t inference with Barnard-Rubin df")
}

charls_rubin_scalar <- function(est, var, label = "quantity", dfcom = Inf,
                                scale = "identity", conf.level = 0.95) {
  if (!is.numeric(est) || !is.numeric(var) || length(est) != length(var)) {
    stop("scalar estimates and variances must be matching numeric vectors", call. = FALSE)
  }
  charls_rubin_pool(lapply(est, function(x) stats::setNames(x, label)),
                    lapply(var, function(x) stats::setNames(x, label)),
                    dfcom = dfcom, scale = scale, conf.level = conf.level)
}
