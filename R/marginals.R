# Marginal distributions for block MINIMA: GEV, Pearson III, Gamma, Weibull.
#
# Orientation, stated once. Minima are negated so that block minima become block
# maxima, which is what extreme-value theory describes; the Jacobian is 1, so
# densities and likelihoods carry over unchanged. Every fit here models discharge
# Q directly: `cdf(q)` is non-exceedance, `ppf(p)` returns a discharge, and the
# T-year low flow is `ppf(1/T)`. The sign flip is confined to the GEV.
#
# Selection is by AIC and the test is a parametric-bootstrap Anderson-Darling.
# A Kolmogorov-Smirnov p-value is reported in a column named KS_p_invalid and
# must not be used to choose between families: with parameters estimated from
# the same sample its rejection rate at alpha = 0.05 is 0.000 and its median
# p-value 0.92, so it cannot discriminate.

suppressPackageStartupMessages({
  library(lmomco)
  library(fitdistrplus)
})

.clip01 <- function(p, eps = 1e-12) pmin(pmax(p, eps), 1 - eps)

# ---------------------------------------------------------------- GEV (minima)

.gev_nllh <- function(par, y) {
  mu <- par[1]; sigma <- par[2]; xi <- par[3]
  if (!is.finite(sigma) || sigma <= 0) return(1e12)
  if (abs(xi) < 1e-8) {
    z <- (y - mu) / sigma
    return(length(y) * log(sigma) + sum(z + exp(-z)))
  }
  t <- 1 + xi * (y - mu) / sigma
  if (any(!is.finite(t)) || any(t <= 0)) return(1e12)
  length(y) * log(sigma) + sum((1 + 1 / xi) * log(t) + t^(-1 / xi))
}

.gev_qx <- function(q, mu, sigma, xi) {
  q <- .clip01(q)
  if (abs(xi) < 1e-8) mu - sigma * log(-log(q))
  else mu + sigma / xi * ((-log(q))^(-xi) - 1)
}

.gev_Fx <- function(x, mu, sigma, xi) {
  if (abs(xi) < 1e-8) return(exp(-exp(-(x - mu) / sigma)))
  t <- 1 + xi * (x - mu) / sigma
  out <- ifelse(t > 0, exp(-t^(-1 / xi)), ifelse(xi > 0, 0, 1))
  as.numeric(out)
}

.gev_fx <- function(x, mu, sigma, xi) {
  if (abs(xi) < 1e-8) {
    z <- (x - mu) / sigma
    return(exp(-z - exp(-z)) / sigma)
  }
  t <- 1 + xi * (x - mu) / sigma
  ifelse(t > 0, t^(-1 / xi - 1) * exp(-t^(-1 / xi)) / sigma, 0)
}

#' GEV fitted to the negated minima, by multistart maximum likelihood.
#'
#' A single-start optimiser with a penalty-based constraint is fragile on block
#' minima: one dominant outlier is enough to leave it at an infeasible parameter
#' vector whose return levels converge to a spurious bound. Several starting
#' values are tried with two optimisers, and a candidate is accepted only if
#' every observation lies inside the fitted support and the objective is finite.
fit_gev_minima <- function(x) {
  x <- as.numeric(x); y <- -x; n <- length(y)
  starts <- list()
  lm <- try(lmomco::pargev(lmomco::lmoms(y)), silent = TRUE)
  if (!inherits(lm, "try-error")) {
    # lmomco reports kappa = -xi in the Coles convention used here
    starts[[length(starts) + 1]] <- c(lm$para[1], lm$para[2], -lm$para[3])
  }
  s <- stats::sd(y); m <- mean(y)
  sig0 <- max(s * sqrt(6) / pi, 1e-6)
  for (xi0 in c(-0.1, 0.1, -0.3, 0.3)) {
    starts[[length(starts) + 1]] <- c(m - 0.5772 * sig0, sig0, xi0)
  }
  best <- NULL
  for (st in starts) {
    for (meth in c("BFGS", "Nelder-Mead")) {
      f <- try(stats::optim(st, .gev_nllh, y = y, method = meth,
                            control = list(maxit = 5000, reltol = 1e-12)), silent = TRUE)
      if (inherits(f, "try-error") || !is.finite(f$value) || f$value >= 1e6) next
      p <- f$par
      t <- 1 + p[3] * (y - p[1]) / p[2]
      if (p[2] <= 0 || any(!is.finite(t)) || any(t <= 0)) next
      if (is.null(best) || f$value < best$value) best <- f
    }
  }
  if (is.null(best)) stop("GEV multistart found no feasible optimum")
  mu <- best$par[1]; sigma <- best$par[2]; xi <- best$par[3]
  notes <- character()
  if (xi < 0) {
    # bounded above in X, i.e. bounded BELOW in Q: physically sensible here
    notes <- c(notes, sprintf("shape %.4f < 0: implied lower bound on Q = %.4g",
                              xi, -(mu - sigma / xi)))
  } else if (xi > 0) {
    notes <- c(notes, sprintf("shape %.4f > 0: Q unbounded below, check the low tail", xi))
  }
  structure(list(
    name = "GEV", n_params = 3L,
    params = c(location = mu, scale = sigma, shape = xi),
    loglik = -best$value, data = x, notes = notes,
    ppf = function(p) -.gev_qx(1 - .clip01(p), mu, sigma, xi),
    cdf = function(q) 1 - .gev_Fx(-q, mu, sigma, xi),
    pdf = function(q) .gev_fx(-q, mu, sigma, xi),
    rand = function(k) -.gev_qx(stats::runif(k), mu, sigma, xi),
    support = c(if (xi < 0) -(mu - sigma / xi) else -Inf, Inf),
    refit = fit_gev_minima
  ), class = "lowflow_marginal")
}

# ------------------------------------------------------------------ Pearson III

#' Pearson type III by L-moments (Hosking), fitted to the discharges directly.
#'
#' Pearson III carries a hard support bound whose position depends on the sign
#' of the fitted skew, and BOTH signs bite at exactly the end of the record a
#' low-flow analysis cares about. Positive skew puts a lower bound on Q: if that
#' bound lands above the sample minimum the density there is zero, the
#' log-likelihood is -Inf and AIC infinite. Negative skew removes the lower
#' bound altogether, and the low tail can then cross zero discharge at a
#' droughtiness level well inside the range people report. Both cases are
#' recorded in `notes`, and the fitted quantile function should always be
#' evaluated at the return period you intend to quote before quoting it.
fit_pe3 <- function(x) {
  x <- as.numeric(x)
  para <- lmomco::parpe3(lmomco::lmoms(x))
  mu <- para$para[1]; sigma <- para$para[2]; gam <- para$para[3]
  lo <- -Inf; hi <- Inf
  if (gam > 0) lo <- mu - 2 * sigma / gam
  if (gam < 0) hi <- mu - 2 * sigma / gam
  notes <- character()
  if (any(x < lo) || any(x > hi)) {
    notes <- c(notes, "sample value outside fitted support: density 0, cdf saturates at 0 or 1")
  }
  if (is.finite(lo) && lo > 0) {
    notes <- c(notes, sprintf("implied lower bound on Q = %.4g", lo))
  }
  if (!is.finite(lo)) {
    q_at_zero <- try(stats::uniroot(function(p) lmomco::quape3(p, para),
                                    c(1e-8, 0.5))$root, silent = TRUE)
    msg <- "skew < 0: no lower bound on Q, the low tail is unbounded"
    if (!inherits(q_at_zero, "try-error")) {
      msg <- sprintf("%s (crosses zero discharge at droughtiness u = %.4f)",
                     msg, 1 - q_at_zero)
    }
    notes <- c(notes, msg)
  }
  dens <- lmomco::pdfpe3(x, para)
  ll <- if (any(!is.finite(dens)) || any(dens <= 0)) -Inf else sum(log(dens))
  structure(list(
    name = "Pearson III", n_params = 3L,
    params = c(mu = mu, sigma = sigma, gamma = gam),
    loglik = ll, data = x, notes = notes,
    ppf = function(p) lmomco::quape3(.clip01(p), para),
    cdf = function(q) lmomco::cdfpe3(q, para),
    pdf = function(q) lmomco::pdfpe3(q, para),
    rand = function(k) lmomco::quape3(stats::runif(k), para),
    support = c(lo, hi), refit = fit_pe3
  ), class = "lowflow_marginal")
}

# ------------------------------------------------------------- Gamma / Weibull

.fit_fitdistr <- function(x, dist, label) {
  x <- as.numeric(x)
  pos <- x[x > 0]
  notes <- character()
  if (length(pos) < length(x)) {
    notes <- c(notes, sprintf("%d zero or negative value(s) excluded from this fit",
                              length(x) - length(pos)))
  }
  f <- fitdistrplus::fitdist(pos, dist)
  e <- f$estimate
  structure(list(
    name = label, n_params = 2L, params = e, loglik = as.numeric(f$loglik),
    data = x, notes = notes,
    ppf = function(p) do.call(paste0("q", dist), c(list(.clip01(p)), as.list(e))),
    cdf = function(q) do.call(paste0("p", dist), c(list(q), as.list(e))),
    pdf = function(q) do.call(paste0("d", dist), c(list(q), as.list(e))),
    rand = function(k) do.call(paste0("r", dist), c(list(k), as.list(e))),
    support = c(0, Inf),
    refit = function(z) .fit_fitdistr(z, dist, label)
  ), class = "lowflow_marginal")
}

fit_gamma <- function(x) .fit_fitdistr(x, "gamma", "Gamma")
fit_weibull <- function(x) .fit_fitdistr(x, "weibull", "Weibull")

#' Fit all four candidate families; failures are reported, never silent.
fit_all <- function(x) {
  fns <- list(GEV = fit_gev_minima, `Pearson III` = fit_pe3,
              Gamma = fit_gamma, Weibull = fit_weibull)
  out <- list()
  for (nm in names(fns)) {
    f <- try(fns[[nm]](x), silent = TRUE)
    out[[nm]] <- if (inherits(f, "try-error")) {
      structure(list(name = nm, failed = TRUE,
                     message = trimws(as.character(f))), class = "lowflow_failed")
    } else f
  }
  out
}

aic_of <- function(fit) 2 * fit$n_params - 2 * fit$loglik

# --------------------------------------------------- goodness of fit, intervals

#' Anderson-Darling A-squared on the probability-integral transform.
ad_statistic <- function(fit, x = NULL) {
  x <- if (is.null(x)) fit$data else x
  z <- sort(.clip01(fit$cdf(sort(x))))
  n <- length(z)
  i <- seq_len(n)
  -n - sum((2 * i - 1) * (log(z) + log(1 - rev(z)))) / n
}

#' Parametric-bootstrap Anderson-Darling test.
#'
#' The null distribution is obtained by simulating from the fitted family and
#' refitting, which is what makes the p-value valid when the parameters were
#' estimated from the same sample.
ad_test <- function(fit, B = 499) {
  obs <- ad_statistic(fit)
  n <- length(fit$data)
  null <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    s <- try(fit$rand(n), silent = TRUE)
    if (inherits(s, "try-error")) next
    g <- try(fit$refit(s), silent = TRUE)
    if (inherits(g, "try-error")) next
    a <- try(ad_statistic(g), silent = TRUE)
    if (!inherits(a, "try-error") && is.finite(a)) null[b] <- a
  }
  null <- null[is.finite(null)]
  p <- if (length(null)) (1 + sum(null >= obs)) / (length(null) + 1) else NA_real_
  list(A2 = obs, p_value = p, B_effective = length(null))
}

#' Kolmogorov-Smirnov p-value with in-sample parameters. Reported, never used.
ks_naive <- function(fit) {
  suppressWarnings(stats::ks.test(fit$data, function(q) fit$cdf(q))$p.value)
}

#' AIC table with the bootstrap Anderson-Darling test and the invalid KS column.
compare_marginals <- function(fits, B = 499) {
  rows <- list()
  for (nm in names(fits)) {
    f <- fits[[nm]]
    if (inherits(f, "lowflow_failed")) {
      rows[[nm]] <- data.frame(distribution = nm, k = NA_integer_, nllh = NA_real_,
                               AIC = NA_real_, dAIC = NA_real_, AD_A2 = NA_real_,
                               AD_p = NA_real_, KS_p_invalid = NA_real_,
                               notes = paste("FIT FAILED:", f$message))
      next
    }
    a <- if (is.finite(f$loglik)) ad_test(f, B) else list(A2 = NA_real_, p_value = NA_real_)
    rows[[nm]] <- data.frame(
      distribution = nm, k = f$n_params, nllh = -f$loglik, AIC = aic_of(f),
      dAIC = NA_real_, AD_A2 = a$A2, AD_p = a$p_value,
      KS_p_invalid = if (is.finite(f$loglik)) ks_naive(f) else NA_real_,
      notes = paste(f$notes, collapse = "; ")
    )
  }
  tab <- do.call(rbind, rows)
  tab$dAIC <- tab$AIC - min(tab$AIC, na.rm = TRUE)
  tab <- tab[order(tab$AIC), ]
  rownames(tab) <- NULL
  tab
}

#' The family with the lowest AIC among those that actually fitted.
select_marginal <- function(fits) {
  ok <- Filter(function(f) !inherits(f, "lowflow_failed") && is.finite(f$loglik), fits)
  if (!length(ok)) stop("no candidate distribution produced a finite likelihood")
  names(ok)[which.min(vapply(ok, aic_of, numeric(1)))]
}

#' T-year low flow, `ppf(1/T)`, for every family that fitted.
return_levels <- function(fits, T = c(2, 5, 10, 20, 50, 100)) {
  out <- data.frame(T = T)
  for (nm in names(fits)) {
    f <- fits[[nm]]
    out[[nm]] <- if (inherits(f, "lowflow_failed")) NA_real_ else f$ppf(1 / T)
  }
  out
}

#' Parametric-bootstrap confidence interval on the T-year low flow.
#'
#' An interval that crosses zero discharge means the fit is extrapolating past
#' physical plausibility, not that negative flow is possible; a warning is
#' emitted when that happens.
return_level_ci <- function(fit, T = c(2, 5, 10, 20, 50, 100), B = 999, level = 0.95) {
  n <- length(fit$data)
  est <- fit$ppf(1 / T)
  boot <- matrix(NA_real_, nrow = B, ncol = length(T))
  for (b in seq_len(B)) {
    s <- try(fit$rand(n), silent = TRUE)
    if (inherits(s, "try-error")) next
    g <- try(fit$refit(s), silent = TRUE)
    if (inherits(g, "try-error")) next
    v <- try(g$ppf(1 / T), silent = TRUE)
    if (!inherits(v, "try-error")) boot[b, ] <- v
  }
  a <- (1 - level) / 2
  lo <- apply(boot, 2, stats::quantile, probs = a, na.rm = TRUE)
  hi <- apply(boot, 2, stats::quantile, probs = 1 - a, na.rm = TRUE)
  if (any(lo < 0, na.rm = TRUE)) {
    warning("a confidence interval crosses zero discharge: the fit is extrapolating ",
            "past physical plausibility. Prefer a lower-bounded family, or fit in log space.")
  }
  data.frame(T = T, distribution = fit$name, estimate = est, lower = lo, upper = hi,
             B_effective = sum(stats::complete.cases(boot)))
}

.gev_ns_nllh <- function(par, y, t) {
  mu0 <- par[1]; mu1 <- par[2]; sigma <- par[3]; xi <- par[4]
  if (!is.finite(sigma) || sigma <= 0) return(1e12)
  mu <- mu0 + mu1 * t
  if (abs(xi) < 1e-8) {
    z <- (y - mu) / sigma
    return(length(y) * log(sigma) + sum(z + exp(-z)))
  }
  w <- 1 + xi * (y - mu) / sigma
  if (any(!is.finite(w)) || any(w <= 0)) return(1e12)
  length(y) * log(sigma) + sum((1 + 1 / xi) * log(w) + w^(-1 / xi))
}

#' Stationary versus trend-in-location GEV, compared by likelihood ratio.
#'
#' A transient projection to 2100 is the case where a stationary GEV is least
#' defensible, so the trend is fitted and tested rather than assumed away. The
#' covariate is the standardised block year, so the reported trend is the shift
#' in the GEV location parameter per standard deviation of year.
#'
#' Both models are fitted by the same multistart search used for the stationary
#' GEV, and the trend model is additionally seeded from the stationary optimum.
#' That seeding matters: the stationary fit is the trend model with mu1 = 0, so
#' the larger model can never have a worse likelihood at a genuine optimum, and
#' a negative deviance is therefore proof of an optimiser failure rather than
#' evidence about the data. Every candidate is also required to keep all
#' observations inside the fitted support.
gev_trend_test <- function(values, years) {
  y <- -as.numeric(values)                      # minima -> maxima
  t <- as.numeric(scale(as.numeric(years)))

  st <- fit_gev_minima(values)                  # stationary, multistart
  p0 <- as.numeric(st$params)
  nllh0 <- -st$loglik

  starts <- list(c(p0[1], 0, p0[2], p0[3]))
  for (m1 in c(-0.5, 0.5, -0.2, 0.2)) {
    starts[[length(starts) + 1]] <- c(p0[1], m1 * p0[2], p0[2], p0[3])
  }
  best <- NULL
  for (s in starts) {
    for (meth in c("BFGS", "Nelder-Mead")) {
      f <- try(stats::optim(s, .gev_ns_nllh, y = y, t = t, method = meth,
                            control = list(maxit = 8000, reltol = 1e-12)), silent = TRUE)
      if (inherits(f, "try-error") || !is.finite(f$value) || f$value >= 1e6) next
      pp <- f$par
      mu <- pp[1] + pp[2] * t
      w <- 1 + pp[4] * (y - mu) / pp[3]
      if (pp[3] <= 0 || any(!is.finite(w)) || any(w <= 0)) next
      if (is.null(best) || f$value < best$value) best <- f
    }
  }
  if (is.null(best)) stop("non-stationary GEV: multistart found no feasible optimum")
  nllh1 <- best$value
  if (nllh1 > nllh0 + 1e-6) {
    stop(sprintf(paste0("non-stationary GEV: the trend model reached a worse likelihood ",
                        "(%.4f) than the nested stationary fit (%.4f), which is an optimiser ",
                        "failure, not a result"), nllh1, nllh0))
  }
  dev <- 2 * (nllh0 - nllh1)
  p <- stats::pchisq(dev, df = 1, lower.tail = FALSE)
  # the trend is fitted on the negated series, so flip the sign back to discharge
  trend_q <- -best$par[2]
  list(
    stationary_params = st$params,
    trend_params = c(mu0 = best$par[1], mu1 = best$par[2],
                     scale = best$par[3], shape = best$par[4]),
    table = data.frame(
      model = c("stationary", "trend in location"),
      k = c(3L, 4L), nllh = c(nllh0, nllh1),
      AIC = c(2 * 3 + 2 * nllh0, 2 * 4 + 2 * nllh1)
    ),
    test = data.frame(
      quantity = c("deviance", "df", "p_value",
                   "location trend per sd of year (discharge units)", "verdict"),
      value = c(sprintf("%.4f", dev), "1", format.pval(p, digits = 4),
                sprintf("%.4f", trend_q),
                if (p < 0.05) "trend supported" else "trend not supported")
    ),
    p_value = p, deviance = dev, location_trend = trend_q
  )
}

#' GEV return levels fitted separately within fixed epochs.
gev_epoch_levels <- function(values, years, epochs, T = c(2, 5, 10, 20, 50, 100)) {
  rows <- list()
  for (e in epochs) {
    sel <- years >= e[1] & years <= e[2]
    if (sum(sel) < 10) next
    f <- try(fit_gev_minima(values[sel]), silent = TRUE)
    if (inherits(f, "try-error")) next
    rows[[length(rows) + 1]] <- data.frame(
      epoch = sprintf("%d-%d", e[1], e[2]), n = sum(sel), T = T, low_flow = f$ppf(1 / T)
    )
  }
  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

#' Trend in the block minima: Mann-Kendall plus a Theil-Sen slope.
trend_tests <- function(values, years) {
  mk <- Kendall::MannKendall(values)
  n <- length(values)
  sl <- numeric(0)
  for (i in seq_len(n - 1)) {
    j <- (i + 1):n
    sl <- c(sl, (values[j] - values[i]) / (years[j] - years[i]))
  }
  ts_slope <- stats::median(sl, na.rm = TRUE)
  lagged <- if (n > 2) stats::acf(values, lag.max = 1, plot = FALSE)$acf[2] else NA_real_
  lb <- if (n > 4) stats::Box.test(values, lag = 1, type = "Ljung-Box")$p.value else NA_real_
  data.frame(
    test = c("Mann-Kendall tau", "Mann-Kendall p", "Theil-Sen slope",
             "lag-1 autocorrelation", "Ljung-Box p (lag 1)"),
    value = c(as.numeric(mk$tau), as.numeric(mk$sl), ts_slope, lagged, lb),
    detail = c("monotone trend in the block minima", "",
               "units per year", "independence between blocks is assumed, not corrected",
               "a significant result means the intervals are too narrow")
  )
}
