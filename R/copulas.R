# Bivariate dependence between two drought series, and joint return periods.
#
# Orientation, stated once. Here `u` is a DROUGHTINESS level, u = 1 - F_Q(q), so
# large u means deep drought. Joint drought is therefore the UPPER set, with
# probability 1 - u - v + C(u, v), and the copulas that represent simultaneous
# drought are those with UPPER-tail dependence. Discharge is recovered from a
# droughtiness level only through discharge_from_u(), which is ppf(1 - u).
#
# Copulas are fitted by maximum pseudo-likelihood on RANK pseudo-observations.
# Ranks are invariant to any increasing marginal transform, so the dependence
# estimate cannot inherit a defect of the fitted margins.
#
# FAMILIES holds only the upper-tail-dependent families, which are the only
# admissible candidates for simultaneous drought: Gaussian, Frank and plain
# Clayton all have lambda_drought = 0, asserting that joint extreme drought
# becomes asymptotically independent. Pass them explicitly for a one-off
# sensitivity check; they do not belong in a candidate set.

suppressPackageStartupMessages(library(copula))

FAMILIES <- c("Gumbel", "Joe", "Survival Clayton")
GAMMA_E <- 0.5772156649015329

#' Rank pseudo-observations, rank / (n + 1), column-wise.
#'
#' Pass the NEGATED minima so that large u means deep drought.
pseudo_obs <- function(...) {
  cols <- lapply(list(...), function(a) as.numeric(a))
  n <- length(cols[[1]])
  if (any(vapply(cols, length, integer(1)) != n)) stop("all series must have the same length")
  do.call(cbind, lapply(cols, function(a) rank(a) / (n + 1)))
}

#' Discharge at a droughtiness level: ppf(1 - u).
discharge_from_u <- function(marg, u) marg$ppf(1 - u)

#' Fit one copula by maximum pseudo-likelihood.
#'
#' Survival Clayton is fitted as a plain Clayton on the reflected ranks 1 - u,
#' which is exactly the survival copula in u-space and avoids any dependence on
#' rotation support in the fitting routine.
fit_copula <- function(family, u) {
  family <- match.arg(family, c("Gumbel", "Joe", "Survival Clayton",
                                "Gaussian", "Frank", "Clayton"))
  n <- nrow(u)
  reflect <- family == "Survival Clayton"
  dat <- if (reflect) 1 - u else u
  base <- switch(family,
    "Gumbel"           = copula::gumbelCopula(2, dim = 2),
    "Joe"              = copula::joeCopula(2, dim = 2),
    "Survival Clayton" = copula::claytonCopula(1, dim = 2),
    "Clayton"          = copula::claytonCopula(1, dim = 2),
    "Frank"            = copula::frankCopula(2, dim = 2),
    "Gaussian"         = copula::normalCopula(0.5, dim = 2)
  )
  f <- copula::fitCopula(base, dat, method = "mpl")
  cop <- f@copula
  theta <- as.numeric(copula::coef(f))[1]
  ll <- as.numeric(stats::logLik(f))

  cdf <- if (reflect) {
    function(a, b) a + b - 1 + copula::pCopula(cbind(1 - a, 1 - b), cop)
  } else {
    function(a, b) copula::pCopula(cbind(a, b), cop)
  }
  dens <- if (reflect) {
    function(a, b) copula::dCopula(cbind(1 - a, 1 - b), cop)
  } else {
    function(a, b) copula::dCopula(cbind(a, b), cop)
  }
  smp <- if (reflect) {
    function(N) 1 - copula::rCopula(N, cop)
  } else {
    function(N) copula::rCopula(N, cop)
  }
  lam <- switch(family,
    "Gumbel"           = 2 - 2^(1 / theta),
    "Joe"              = 2 - 2^(1 / theta),
    "Survival Clayton" = 2^(-1 / theta),
    0
  )
  structure(list(family = family, theta = theta, loglik = ll, n = n,
                 aic = 2 * 1 - 2 * ll, tau = copula::tau(cop),
                 lambda_drought = lam, cop = cop, reflected = reflect,
                 data = u, cdf = cdf, dens = dens, sample = smp),
            class = "lowflow_copula")
}

fit_all_copulas <- function(u, families = FAMILIES) {
  out <- list()
  for (fam in families) {
    f <- try(fit_copula(fam, u), silent = TRUE)
    out[[fam]] <- if (inherits(f, "try-error")) {
      structure(list(family = fam, failed = TRUE,
                     message = trimws(as.character(f))), class = "lowflow_failed")
    } else f
  }
  out
}

#' Parametric-bootstrap Cramer-von Mises (Sn) goodness of fit.
#'
#' B = 499 is the practical minimum: the Monte-Carlo standard error on a p-value
#' near 0.5 at B = 100 is about +/-0.05, too coarse to separate candidates.
gof_sn <- function(model, B = 499) {
  dat <- if (model$reflected) 1 - model$data else model$data
  g <- try(copula::gofCopula(model$cop, dat, N = B, method = "Sn",
                             simulation = "pb", estim.method = "mpl",
                             verbose = FALSE), silent = TRUE)
  if (inherits(g, "try-error")) {
    return(list(copula = model$family, Sn = NA_real_, p_value = NA_real_, B = B))
  }
  list(copula = model$family, Sn = as.numeric(g$statistic),
       p_value = as.numeric(g$p.value), B = B)
}

#' Rank the candidates by AIC and test each one.
compare_copulas <- function(fits, B = 499) {
  rows <- list()
  for (nm in names(fits)) {
    f <- fits[[nm]]
    if (inherits(f, "lowflow_failed")) {
      rows[[nm]] <- data.frame(copula = nm, theta = NA_real_, tau = NA_real_,
                               lambda_drought = NA_real_, loglik = NA_real_,
                               AIC = NA_real_, dAIC = NA_real_, Sn = NA_real_,
                               Sn_p = NA_real_,
                               notes = paste("FIT FAILED:", f$message))
      next
    }
    g <- gof_sn(f, B)
    rows[[nm]] <- data.frame(copula = nm, theta = f$theta, tau = f$tau,
                             lambda_drought = f$lambda_drought, loglik = f$loglik,
                             AIC = f$aic, dAIC = NA_real_, Sn = g$Sn,
                             Sn_p = g$p_value, notes = "")
  }
  tab <- do.call(rbind, rows)
  tab$dAIC <- tab$AIC - min(tab$AIC, na.rm = TRUE)
  tab <- tab[order(tab$AIC), ]
  rownames(tab) <- NULL
  tab
}

select_copula <- function(fits) {
  ok <- Filter(function(f) !inherits(f, "lowflow_failed"), fits)
  if (!length(ok)) stop("no copula could be fitted")
  names(ok)[which.min(vapply(ok, function(f) f$aic, numeric(1)))]
}

#' Caperaa-Fougeres-Genest Pickands estimate; lambda_drought = 2 - 2 A(1/2).
pickands_cfg <- function(u, t = 0.5) {
  a <- -log(pmin(pmax(u[, 1], 1e-12), 1 - 1e-12))
  b <- -log(pmin(pmax(u[, 2], 1e-12), 1 - 1e-12))
  xi <- pmin(a / (1 - t), b / t)
  exp(-GAMMA_E - mean(log(xi)))
}

lambda_cfg <- function(u) 2 - 2 * pickands_cfg(u, 0.5)

#' Model-free dependence measures, so the copula choice has something to answer to.
dependence_summary <- function(u) {
  a <- u[, 1]; b <- u[, 2]
  kt <- suppressWarnings(stats::cor.test(a, b, method = "kendall"))
  sp <- suppressWarnings(stats::cor.test(a, b, method = "spearman"))
  rows <- list(
    data.frame(measure = "Kendall tau", value = as.numeric(kt$estimate),
               p_value = kt$p.value,
               detail = "rank correlation of the two block-minimum series"),
    data.frame(measure = "Spearman rho", value = as.numeric(sp$estimate),
               p_value = sp$p.value, detail = "")
  )
  for (k in c(0.1, 0.2, 0.3)) {
    thr <- 1 - k
    sel <- a > thr
    rows[[length(rows) + 1]] <- data.frame(
      measure = sprintf("empirical lambda_drought (k=%.1f)", k),
      value = if (any(sel)) mean(b[sel] > thr) else NA_real_,
      p_value = NA_real_,
      detail = sprintf("P(both in driest %d%% | one is), n=%d", round(k * 100), sum(sel))
    )
  }
  rows[[length(rows) + 1]] <- data.frame(
    measure = "lambda_drought (CFG)", value = lambda_cfg(u), p_value = NA_real_,
    detail = "non-parametric; compare against the fitted copula's lambda_drought"
  )
  do.call(rbind, rows)
}

# -------------------------------------------------------- joint return periods

#' P(both rivers below their thresholds).
and_probability <- function(model, u, v) 1 - u - v + model$cdf(u, v)

#' P(at least one river below its threshold).
or_probability <- function(model, u, v) 1 - model$cdf(u, v)

.solve_u <- function(f, lo = 1e-9, hi = 1 - 1e-9) {
  r <- try(stats::uniroot(f, c(lo, hi), tol = 1e-12)$root, silent = TRUE)
  if (inherits(r, "try-error")) NA_real_ else r
}

#' The equal-droughtiness point (u = v) on the AND curve of return period T.
#'
#' One point on the curve, not the only one; see most_likely_design_point().
and_design_point <- function(model, T, marg_x, marg_y) {
  rows <- lapply(T, function(Ti) {
    w <- .solve_u(function(w) and_probability(model, w, w) - 1 / Ti)
    data.frame(T = Ti, u = w,
               x = if (is.na(w)) NA_real_ else discharge_from_u(marg_x, w),
               y = if (is.na(w)) NA_real_ else discharge_from_u(marg_y, w))
  })
  do.call(rbind, rows)
}

#' The equal-droughtiness point on the OR contour of return period T.
or_design_point <- function(model, T, marg_x, marg_y) {
  rows <- lapply(T, function(Ti) {
    w <- .solve_u(function(w) or_probability(model, w, w) - 1 / Ti)
    data.frame(T = Ti, u = w,
               x = if (is.na(w)) NA_real_ else discharge_from_u(marg_x, w),
               y = if (is.na(w)) NA_real_ else discharge_from_u(marg_y, w))
  })
  do.call(rbind, rows)
}

#' The AND contour, solved by root-finding rather than filtered off a grid.
#'
#' For each u the matching v is found with uniroot, so every point lies exactly
#' on the contour and the curve comes out monotone and evenly sampled.
and_curve <- function(model, T, marg_x, marg_y, n = 200) {
  target <- 1 / T
  eps <- 1e-9
  # The contour exists only where the AND probability can reach the target at
  # all: P(U>u, V>v) <= 1 - u, so u must be at or below 1 - target. Scanning the
  # whole unit interval and keeping the u values that bracket a root covers that
  # automatically, and never silently returns an empty curve for small T.
  us <- seq(eps, 1 - eps, length.out = n)
  rows <- lapply(us, function(uu) {
    f <- function(v) and_probability(model, uu, v) - target
    if (!is.finite(f(eps)) || !is.finite(f(1 - eps))) return(NULL)
    if (f(eps) * f(1 - eps) > 0) return(NULL)
    vv <- .solve_u(f, eps, 1 - eps)
    if (is.na(vv)) return(NULL)
    data.frame(T = T, u = uu, v = vv,
               x = discharge_from_u(marg_x, uu), y = discharge_from_u(marg_y, vv))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) return(data.frame())
  out[order(out$u), ]
}

#' The most likely realisation on an AND contour.
#'
#' Maximises the joint density along the contour, which answers "if this joint
#' drought happens, what does it most probably look like" -- a different and
#' usually more useful question than the equal-droughtiness point.
most_likely_design_point <- function(model, T, marg_x, marg_y, n = 400) {
  cv <- and_curve(model, T, marg_x, marg_y, n = n)
  if (!nrow(cv)) return(data.frame(T = T, u = NA_real_, v = NA_real_,
                                   x = NA_real_, y = NA_real_))
  fx <- marg_x$pdf(cv$x); fy <- marg_y$pdf(cv$y)
  cd <- model$dens(cv$u, cv$v)
  dens <- cd * fx * fy
  i <- which.max(ifelse(is.finite(dens), dens, -Inf))
  data.frame(T = T, u = cv$u[i], v = cv$v[i], x = cv$x[i], y = cv$y[i])
}

#' AND-probability level whose critical layer has Kendall return period T.
#'
#' Monte Carlo from the fitted copula, so these carry their own sampling error,
#' unlike the analytic AND and OR quantities.
kendall_critical_level <- function(model, T, N = 40000) {
  s <- model$sample(N)
  w <- and_probability(model, s[, 1], s[, 2])
  rows <- lapply(T, function(Ti) {
    q <- as.numeric(stats::quantile(w, 1 / Ti, na.rm = TRUE))
    data.frame(T = Ti, and_probability_level = q,
               T_and_equivalent = if (q > 0) 1 / q else Inf)
  })
  do.call(rbind, rows)
}

#' Return period of the critical layer at AND-probability level p.
kendall_return_period <- function(model, p, N = 40000) {
  s <- model$sample(N)
  w <- and_probability(model, s[, 1], s[, 2])
  pr <- mean(w <= p, na.rm = TRUE)
  if (pr == 0) Inf else 1 / pr
}
