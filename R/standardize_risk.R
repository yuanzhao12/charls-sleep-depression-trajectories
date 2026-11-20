# Conditional-on-class standardization over a fixed empirical target population.
charls_standardize_risk <- function(fit, data, exposure, exposure_levels,
                                    cluster = NULL) {
  if (!inherits(fit, "glm") || fit$family$family != "binomial" || fit$family$link != "logit") {
    stop("standardization requires a logistic glm", call. = FALSE)
  }
  if (missing(exposure_levels) || length(exposure_levels) != 2L ||
      anyNA(exposure_levels) || !exposure %in% names(data) || nrow(data) == 0) {
    stop("explicit reference/exposed levels and nonempty data required", call. = FALSE)
  }
  beta <- stats::coef(fit)
  if (any(!is.finite(beta))) stop("nonfinite or aliased coefficients", call. = FALSE)
  if (length(attr(stats::terms(fit), "offset"))) stop("offset models are not supported", call. = FALSE)
  lev <- fit$xlevels[[exposure]]
  if (is.null(lev)) lev <- levels(data[[exposure]])
  if (is.null(lev) || any(!exposure_levels %in% lev)) stop("unknown exposure level", call. = FALSE)
  terms <- stats::delete.response(stats::terms(fit))
  scenario <- function(level) {
    d <- data
    d[[exposure]] <- factor(rep(level, nrow(d)), levels = lev)
    mf <- stats::model.frame(terms, d, na.action = stats::na.fail, xlev = fit$xlevels)
    X <- stats::model.matrix(terms, mf, contrasts.arg = fit$contrasts)
    if (!setequal(colnames(X), names(beta)) || nrow(X) != nrow(data)) {
      stop("scenario design matrix does not match model", call. = FALSE)
    }
    X <- X[, names(beta), drop = FALSE]
    p <- as.numeric(stats::plogis(X %*% beta))
    list(risk = mean(p), gradient = colMeans(X * (p * (1 - p))), X = X)
  }
  ref <- scenario(exposure_levels[1]); alt <- scenario(exposure_levels[2])
  if (ref$risk <= 0 || alt$risk <= 0) stop("log risk ratio undefined at zero risk", call. = FALSE)
  G <- rbind(risk_ref = ref$gradient, risk_exp = alt$gradient,
             rd = alt$gradient - ref$gradient,
             log_rr = alt$gradient / alt$risk - ref$gradient / ref$risk)
  V <- if (is.null(cluster)) stats::vcov(fit) else charls_cluster_vcov(fit, cluster)$vcov
  V <- charls_check_vcov(V, names(beta))
  U <- G %*% V %*% t(G)
  if (any(diag(U) < -1e-10)) stop("negative derived variance", call. = FALSE)
  estimates <- c(risk_ref = ref$risk, risk_exp = alt$risk,
                 rd = alt$risk - ref$risk, log_rr = log(alt$risk / ref$risk))
  se <- stats::setNames(sqrt(pmax(diag(U), 0)), rownames(G))
  list(risk_ref = ref$risk, risk_exp = alt$risk, rd = unname(estimates["rd"]),
       rr = exp(unname(estimates["log_rr"])), log_rr = unname(estimates["log_rr"]),
       estimates = estimates, se = se, gradient = G, vcov = U, vcov_beta = V,
       n = nrow(data), exposure_levels = exposure_levels,
       target = "fixed empirical covariate distribution of the same analysis frame",
       conditional_on = c("fitted trajectory classification", "empirical target covariates"),
       uncertainty_includes = "model coefficient covariance; MI variability added during pooling",
       uncertainty_excludes = "refitting trajectory model; resampling target distribution")
}

