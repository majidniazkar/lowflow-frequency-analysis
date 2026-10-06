# One univariate and one bivariate implementation, driven by a configuration
# list. Variants are produced by passing a different configuration, never by
# copying the pipeline: the block definition, the output directory and the
# figure labels then exist in exactly one place each and cannot drift apart.

.logger <- function(path) {
  con <- file(path, open = "wt")
  list(
    msg = function(...) {
      txt <- paste0(...)
      cat(txt, "\n", sep = "")
      writeLines(txt, con)
    },
    table = function(df) {
      txt <- utils::capture.output(print(df, row.names = FALSE, digits = 6))
      cat(txt, sep = "\n"); cat("\n")
      writeLines(txt, con)
    },
    close = function() close(con)
  )
}

.csv <- function(df, outdir, stem, written) {
  p <- file.path(outdir, paste0(stem, ".csv"))
  utils::write.csv(df, p, row.names = FALSE)
  c(written, p)
}

univariate_config <- function(path, station, outdir, label = station, season = NULL,
                              D = 7, min_coverage = 0.9, start_month = 1L,
                              zero_policy = "keep", unit = "m^3/s",
                              T = c(2, 5, 10, 20, 50, 100), B_gof = 499, B_ci = 999,
                              nonstationary = FALSE, epochs = list(),
                              scan_season = TRUE, date_col = "date", flow_col = "flow") {
  as.list(environment())
}

bivariate_config <- function(path_x, path_y, name_x, name_y, outdir,
                             label = paste(name_x, "and", name_y), season = NULL,
                             D = 7, min_coverage = 0.9, start_month = 1L,
                             zero_policy = "keep", unit = "m^3/s",
                             T = c(2, 5, 10, 20, 50, 100), B_gof = 499,
                             families = FAMILIES, N_kendall = 40000,
                             date_col = "date", flow_col = "flow") {
  as.list(environment())
}

#' Full marginal low-flow analysis for one station.
run_univariate <- function(cfg) {
  dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
  log <- .logger(file.path(cfg$outdir, "run_log.txt"))
  on.exit(log$close(), add = TRUE)
  written <- character()

  log$msg(sprintf("[1] Record: %s", cfg$label))
  rec <- read_flow(cfg$path, cfg$station, cfg$D, cfg$zero_policy, cfg$date_col, cfg$flow_col)
  out <- utils::capture.output(summary(rec))
  log$msg(paste(out, collapse = "\n"))

  season <- cfg$season
  if (cfg$scan_season) {
    scan <- scan_seasons(rec, min_coverage = cfg$min_coverage, start_month = cfg$start_month)
    log$msg("\n[2] Candidate low-flow windows, ranked by capture of the annual minimum")
    log$table(utils::head(scan[, setdiff(names(scan), "months")], 10))
    written <- .csv(scan, cfg$outdir, "season_scan", written)
    lab <- if (is.null(season)) NULL else season_label(season)
    if (!is.null(season)) {
      row <- scan[scan$window == lab, ]
      if (nrow(row)) {
        log$msg(sprintf("    chosen window %s: capture rate %.1f%%", lab, 100 * row$capture_rate[1]))
        if (row$capture_rate[1] < 0.8) {
          log$msg(sprintf("    WARNING: %s contains the annual minimum in only %.0f%% of blocks.",
                          lab, 100 * row$capture_rate[1]))
        }
      }
    }
    written <- c(written, save_fig(
      fig_season_scan(scan, lab, sprintf("%s: does the window contain the drought?", cfg$station)),
      cfg$outdir, "fig_season_scan", height = 5.2))
  }

  bm <- block_minima(rec, season, cfg$min_coverage, cfg$start_month)
  x <- bm$values
  log$msg(sprintf("\n[3] Block minima: %d used, %d rejected on coverage (min_coverage = %.2f)",
                  nrow(bm$used), nrow(bm$rejected), cfg$min_coverage))
  if (nrow(bm$rejected)) log$table(bm$rejected[, c("year", "coverage", "n_valid", "n_possible")])
  if (length(x) < 15) log$msg("    WARNING: fewer than 15 blocks; intervals will be very wide.")
  written <- .csv(bm$table, cfg$outdir, "block_minima", written)

  tr <- trend_tests(x, bm$years)
  log$msg("\n[4] Trend and independence of the block minima")
  log$table(tr)
  written <- .csv(tr, cfg$outdir, "trend_tests", written)
  slope <- tr$value[tr$test == "Theil-Sen slope"]
  mkp <- tr$value[tr$test == "Mann-Kendall p"]
  written <- c(written, save_fig(
    fig_block_minima(bm, sprintf("%s: block minima of the %d-day mean", cfg$station, cfg$D),
                     cfg$unit, slope, mkp), cfg$outdir, "fig_block_minima"))

  log$msg("\n[5] Candidate distributions, ranked by AIC")
  fits <- fit_all(x)
  cmp <- compare_marginals(fits, cfg$B_gof)
  log$table(cmp)
  log$msg("    KS_p_invalid is printed to show why it must not be used for selection:")
  log$msg("    with in-sample parameters its rejection rate at alpha = 0.05 is 0.000.")
  written <- .csv(cmp, cfg$outdir, "distribution_comparison", written)
  best <- select_marginal(fits)
  log$msg(sprintf("    selected by AIC: %s", best))
  for (nm in names(fits)) {
    f <- fits[[nm]]
    if (!inherits(f, "lowflow_failed") && length(f$notes)) {
      log$msg(sprintf("    note (%s): %s", nm, paste(f$notes, collapse = "; ")))
    }
  }
  written <- c(written, save_fig(fig_density_fits(
    x, fits, sprintf("%s: fitted densities", cfg$label), cfg$unit),
    cfg$outdir, "fig_density_fits"))
  written <- c(written, save_fig(fig_qq(
    x, fits, sprintf("%s: Q-Q", cfg$label), cfg$unit),
    cfg$outdir, "fig_qq", height = 5.4))

  log$msg("\n[6] Low-flow return levels")
  rl <- return_levels(fits, cfg$T)
  log$table(rl)
  written <- .csv(rl, cfg$outdir, "return_levels_all_distributions", written)
  ci <- return_level_ci(fits[[best]], cfg$T, cfg$B_ci)
  log$msg(sprintf("    95%% intervals on %s (parametric bootstrap, B = %d)", best, cfg$B_ci))
  log$table(ci)
  written <- .csv(ci, cfg$outdir, "return_levels_ci_selected", written)
  written <- c(written, save_fig(fig_return_levels(
    x, fits, best, ci, sprintf("%s: low-flow return levels", cfg$label), cfg$unit),
    cfg$outdir, "fig_return_levels"))

  ns <- NULL
  if (isTRUE(cfg$nonstationary)) {
    log$msg("\n[7] Stationary versus trend-in-location GEV")
    ns <- try(gev_trend_test(x, bm$years), silent = TRUE)
    if (inherits(ns, "try-error")) {
      log$msg(sprintf("    non-stationary fit failed: %s", trimws(as.character(ns))))
      ns <- NULL
    } else {
      log$table(ns$table)
      log$table(ns$test)
      written <- .csv(ns$table, cfg$outdir, "nonstationary_fit", written)
      written <- .csv(ns$test, cfg$outdir, "nonstationary_test", written)
    }
    if (length(cfg$epochs)) {
      ep <- gev_epoch_levels(x, bm$years, cfg$epochs, cfg$T)
      if (nrow(ep)) {
        log$msg("    GEV return levels fitted separately by epoch")
        log$table(ep)
        written <- .csv(ep, cfg$outdir, "return_levels_by_epoch", written)
      }
    }
    log$msg(paste0("    NOTE: a projected low flow is reportable as a RELATIVE change only, ",
                   "until the simulated\n    series has been bias-adjusted against the observed ",
                   "record in its LOW TAIL, and until it\n    is stated whether the file is one ",
                   "realisation or one member of an ensemble."))
  }

  log$msg(sprintf("\n%d files written to %s", length(written), cfg$outdir))
  invisible(list(record = rec, blocks = bm, fits = fits, comparison = cmp,
                 selected = best, return_levels = rl, ci = ci, trend = tr,
                 nonstationary = ns, season = season, written = written,
                 outdir = cfg$outdir))
}

#' Pair two records on days valid at BOTH gauges, over their full common span.
#'
#' Non-common days are masked to NA rather than dropped, so the coverage screen
#' still divides by the number of calendar days the block spans. Reindexing onto
#' the common days alone would make every coverage exactly 1 and disable the
#' screen without any error.
joint_records <- function(rec_x, rec_y) {
  dx <- rec_x$daily; dy <- rec_y$daily
  full <- seq(min(c(dx$date, dy$date)), max(c(dx$date, dy$date)), by = "day")
  ax <- rep(NA_real_, length(full)); ay <- ax
  ax[match(dx$date, full)] <- dx$flow_D
  ay[match(dy$date, full)] <- dy$flow_D
  common <- !is.na(ax) & !is.na(ay)
  ax[!common] <- NA_real_
  ay[!common] <- NA_real_
  mk <- function(a, rec) {
    r <- rec
    r$daily <- data.frame(date = full, flow = a, flow_D = a)
    r
  }
  list(x = mk(ax, rec_x), y = mk(ay, rec_y), n_common = sum(common),
       span = range(full))
}

#' Joint low-flow analysis for two gauges.
run_bivariate <- function(cfg) {
  dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
  log <- .logger(file.path(cfg$outdir, "run_log.txt"))
  on.exit(log$close(), add = TRUE)
  written <- character()

  log$msg(sprintf("[1] Records: %s", cfg$label))
  ra <- read_flow(cfg$path_x, cfg$name_x, cfg$D, cfg$zero_policy, cfg$date_col, cfg$flow_col)
  rb <- read_flow(cfg$path_y, cfg$name_y, cfg$D, cfg$zero_policy, cfg$date_col, cfg$flow_col)
  log$msg(paste(utils::capture.output(summary(ra)), collapse = "\n"))
  log$msg(paste(utils::capture.output(summary(rb)), collapse = "\n"))

  jr <- joint_records(ra, rb)
  log$msg(sprintf("\n[2] Days with a valid %d-day mean at BOTH gauges: %d (%s to %s)",
                  cfg$D, jr$n_common, jr$span[1], jr$span[2]))
  if (jr$n_common == 0) {
    stop("No day has a valid D-day mean at both gauges, so no block can be paired. ",
         "Check that the two records overlap in time and that both parse to dates.")
  }

  bma <- block_minima(jr$x, cfg$season, cfg$min_coverage, cfg$start_month)
  bmb <- block_minima(jr$y, cfg$season, cfg$min_coverage, cfg$start_month)
  yrs <- intersect(bma$years, bmb$years)
  if (!length(yrs)) {
    stop("No block is usable for BOTH gauges after the coverage screen. ",
         sprintf("%s kept %d blocks, %s kept %d, with no year in common.",
                 cfg$name_x, nrow(bma$used), cfg$name_y, nrow(bmb$used)))
  }
  qx <- bma$used$value[match(yrs, bma$years)]
  qy <- bmb$used$value[match(yrs, bmb$years)]
  log$msg(sprintf("    paired blocks: %d (%d-%d)", length(yrs), min(yrs), max(yrs)))
  if (length(yrs) < 20) {
    log$msg(sprintf("    WARNING: only %d paired blocks. Do not report beyond T = %d.",
                    length(yrs), max(2, round(2 * length(yrs) / 10) * 10)))
  }
  paired <- data.frame(year = yrs, x = qx, y = qy)
  names(paired) <- c("year", cfg$name_x, cfg$name_y)
  written <- .csv(paired, cfg$outdir, "paired_block_minima", written)

  log$msg("\n[3] Marginals, each gauge on its own evidence")
  fa <- fit_all(qx); fb <- fit_all(qy)
  ca <- compare_marginals(fa, cfg$B_gof); cb <- compare_marginals(fb, cfg$B_gof)
  log$msg(sprintf("  %s:", cfg$name_x)); log$table(ca)
  log$msg(sprintf("  %s:", cfg$name_y)); log$table(cb)
  ba <- select_marginal(fa); bb <- select_marginal(fb)
  log$msg(sprintf("    selected: %s -> %s, %s -> %s", cfg$name_x, ba, cfg$name_y, bb))
  written <- .csv(ca, cfg$outdir, paste0("marginal_comparison_", cfg$name_x), written)
  written <- .csv(cb, cfg$outdir, paste0("marginal_comparison_", cfg$name_y), written)
  mx <- fa[[ba]]; my <- fb[[bb]]

  log$msg("\n[4] Dependence of the two drought series (model-free)")
  u <- pseudo_obs(-qx, -qy)
  dep <- dependence_summary(u)
  log$table(dep)
  written <- .csv(dep, cfg$outdir, "dependence_summary", written)
  tau_emp <- dep$value[dep$measure == "Kendall tau"]
  lam_emp <- dep$value[dep$measure == "lambda_drought (CFG)"]
  written <- c(written, save_fig(fig_pseudo_obs(
    u, cfg$name_x, cfg$name_y, tau_emp, lam_emp,
    sprintf("%s: droughtiness pseudo-observations", cfg$label)),
    cfg$outdir, "fig_pseudo_obs", width = 6, height = 6))

  log$msg("\n[5] Copulas, fitted on ranks by maximum pseudo-likelihood, ranked by AIC")
  cops <- fit_all_copulas(u, cfg$families)
  ccmp <- compare_copulas(cops, cfg$B_gof)
  log$table(ccmp)
  written <- .csv(ccmp, cfg$outdir, "copula_comparison", written)
  bestc <- select_copula(cops)
  model <- cops[[bestc]]
  log$msg(sprintf("    selected by AIC: %s (theta = %.6f, tau = %.4f, lambda_drought = %.4f)",
                  bestc, model$theta, model$tau, model$lambda_drought))
  log$msg(sprintf("    empirical tau = %.4f, CFG lambda_drought = %.4f", tau_emp, lam_emp))
  log$msg(paste0("    All candidates here have upper-tail dependence on the droughtiness\n",
                 "    scale. Report the lambda_drought range across them rather than ",
                 "defending one AIC win."))

  log$msg("\n[6] Joint return levels. AND, OR and Kendall are three different events.")
  dp <- and_design_point(model, cfg$T, mx, my)
  names(dp) <- c("T", "u", cfg$name_x, cfg$name_y)
  log$msg("  AND: both gauges below their thresholds, at the equal-droughtiness point")
  log$table(dp)
  written <- .csv(dp, cfg$outdir, "and_design_points", written)

  op <- or_design_point(model, cfg$T, mx, my)
  names(op) <- c("T", "u", cfg$name_x, cfg$name_y)
  log$msg("  OR: at least one gauge below its threshold")
  log$table(op)
  written <- .csv(op, cfg$outdir, "or_design_points", written)

  mld <- do.call(rbind, lapply(cfg$T, function(t) most_likely_design_point(model, t, mx, my)))
  log$msg("  Most likely realisation on each AND contour")
  log$table(mld)
  written <- .csv(mld, cfg$outdir, "and_most_likely_realisations", written)

  kc <- kendall_critical_level(model, cfg$T, cfg$N_kendall)
  log$msg(sprintf("  Kendall critical layers (Monte Carlo, N = %d, carries sampling error)",
                  cfg$N_kendall))
  log$table(kc)
  written <- .csv(kc, cfg$outdir, "kendall_critical_levels", written)

  curves <- do.call(rbind, lapply(cfg$T, function(t) and_curve(model, t, mx, my, n = 240)))
  written <- .csv(curves, cfg$outdir, "and_curves", written)
  obs <- data.frame(x = qx, y = qy)
  written <- c(written, save_fig(fig_and_curves(
    curves, dp[, c("T", "u", cfg$name_x, cfg$name_y)] |>
      setNames(c("T", "u", "x", "y")), mld, obs,
    cfg$name_x, cfg$name_y, cfg$label, cfg$unit, bestc),
    cfg$outdir, "fig_and_curves", width = 7.8, height = 5))

  runner <- setdiff(ccmp$copula[!is.na(ccmp$AIC)], bestc)
  if (length(runner)) {
    m2 <- cops[[runner[1]]]
    cv2 <- do.call(rbind, lapply(cfg$T, function(t) and_curve(m2, t, mx, my, n = 240)))
    dp2 <- and_design_point(m2, cfg$T, mx, my)
    written <- c(written, save_fig(fig_and_curves(
      cv2, dp2, NULL, obs, cfg$name_x, cfg$name_y,
      paste0(cfg$label, " [runner-up by AIC, shown for sensitivity]"),
      cfg$unit, runner[1]), cfg$outdir, "fig_and_curves_runner_up",
      width = 7.8, height = 5))
    sens <- do.call(rbind, lapply(names(Filter(function(f) !inherits(f, "lowflow_failed"), cops)),
      function(nm) {
        d <- and_design_point(cops[[nm]], cfg$T, mx, my)
        data.frame(copula = nm, T = d$T, u = d$u, x = d$x, y = d$y)
      }))
    names(sens) <- c("copula", "T", "u", cfg$name_x, cfg$name_y)
    log$msg("  Sensitivity of the AND design point to the copula choice")
    log$table(sens)
    written <- .csv(sens, cfg$outdir, "and_design_point_copula_sensitivity", written)
  }

  log$msg(sprintf("\n%d files written to %s", length(written), cfg$outdir))
  invisible(list(paired = paired, marginals = list(x = fa, y = fb),
                 selected_marginals = c(ba, bb), u = u, dependence = dep,
                 copulas = cops, copula_comparison = ccmp, selected_copula = bestc,
                 and_design = dp, or_design = op, most_likely = mld,
                 kendall = kc, curves = curves, written = written, outdir = cfg$outdir))
}
