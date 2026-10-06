# Figures. Every one is generated from the same fitted objects the tables are
# built from, so a label can never describe a different fit than the curve it
# sits on, and a seasonal curve can never be drawn against annual points.

suppressPackageStartupMessages(library(ggplot2))

PALETTE <- c("#482878", "#1F9E63", "#D1495B", "#2A7FA8", "#E0A526", "#6E6E6E")

theme_lowflow <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey90", linewidth = 0.3),
      plot.title = ggplot2::element_text(size = base_size * 1.05, face = "plain"),
      plot.subtitle = ggplot2::element_text(size = base_size * 0.9, colour = "grey30"),
      legend.position = "right", legend.title = ggplot2::element_blank(),
      axis.title = ggplot2::element_text(size = base_size)
    )
}

#' Block minima through time, with rejected blocks drawn as open markers.
fig_block_minima <- function(bm, title = "", unit = "m^3/s", slope = NA, mk_p = NA) {
  used <- bm$used; rej <- bm$rejected
  sub <- if (is.finite(slope)) {
    sprintf("Theil-Sen %.4g %s/yr (Mann-Kendall p = %.3f)", slope, unit, mk_p)
  } else ""
  p <- ggplot2::ggplot(used, ggplot2::aes(year, value)) +
    ggplot2::geom_line(colour = "grey20", linewidth = 0.4) +
    ggplot2::geom_point(colour = "black", size = 1.6)
  if (nrow(rej)) {
    p <- p + ggplot2::geom_point(data = rej, shape = 21, fill = NA,
                                 colour = PALETTE[3], size = 2.2, stroke = 0.8)
  }
  if (is.finite(slope)) {
    p <- p + ggplot2::geom_abline(
      intercept = mean(used$value) - slope * mean(used$year),
      slope = slope, colour = PALETTE[3], linetype = "dashed", linewidth = 0.6)
  }
  p + ggplot2::labs(
    title = title, subtitle = sub, x = "Year",
    y = sprintf("Minimum discharge (%s)", unit),
    caption = sprintf("filled = used (n = %d); open = rejected on coverage (n = %d)",
                      nrow(used), nrow(rej))) +
    theme_lowflow()
}

#' Fitted densities against the histogram of the block minima.
fig_density_fits <- function(x, fits, title = "", unit = "m^3/s") {
  ok <- Filter(function(f) !inherits(f, "lowflow_failed"), fits)
  grid <- seq(max(0, min(x) * 0.5), max(x) * 1.25, length.out = 400)
  dd <- do.call(rbind, lapply(names(ok), function(nm) {
    d <- try(ok[[nm]]$pdf(grid), silent = TRUE)
    if (inherits(d, "try-error")) return(NULL)
    data.frame(q = grid, density = as.numeric(d), distribution = nm)
  }))
  ggplot2::ggplot() +
    ggplot2::geom_histogram(data = data.frame(x = x), ggplot2::aes(x, after_stat(density)),
                            bins = 10, fill = "grey85", colour = "white") +
    ggplot2::geom_line(data = dd, ggplot2::aes(q, density, colour = distribution),
                       linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = PALETTE) +
    ggplot2::labs(title = title, x = sprintf("Block-minimum discharge (%s)", unit),
                  y = "Probability density",
                  caption = sprintf("observed n = %d", length(x))) +
    theme_lowflow()
}

#' Q-Q plot, one panel per candidate.
fig_qq <- function(x, fits, title = "", unit = "m^3/s") {
  ok <- Filter(function(f) !inherits(f, "lowflow_failed"), fits)
  n <- length(x)
  pp <- (seq_len(n) - 0.44) / (n + 0.12)   # Gringorten plotting positions
  xs <- sort(x)
  dd <- do.call(rbind, lapply(names(ok), function(nm) {
    th <- try(ok[[nm]]$ppf(pp), silent = TRUE)
    if (inherits(th, "try-error")) return(NULL)
    data.frame(theoretical = as.numeric(th), observed = xs, distribution = nm)
  }))
  ggplot2::ggplot(dd, ggplot2::aes(theoretical, observed)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey50") +
    ggplot2::geom_point(ggplot2::aes(colour = distribution), size = 1.4, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = PALETTE) +
    ggplot2::facet_wrap(~ distribution) +
    ggplot2::labs(title = title, x = sprintf("Theoretical quantile (%s)", unit),
                  y = sprintf("Observed quantile (%s)", unit)) +
    theme_lowflow()
}

#' Low-flow return levels with a confidence band on the selected family.
fig_return_levels <- function(x, fits, selected, ci, title = "", unit = "m^3/s") {
  ok <- Filter(function(f) !inherits(f, "lowflow_failed"), fits)
  Tg <- exp(seq(log(1.5), log(max(ci$T) * 1.2), length.out = 60))
  curves <- do.call(rbind, lapply(names(ok), function(nm) {
    v <- try(ok[[nm]]$ppf(1 / Tg), silent = TRUE)
    if (inherits(v, "try-error")) return(NULL)
    data.frame(T = Tg, q = as.numeric(v), distribution = nm)
  }))
  n <- length(x)
  pp <- (seq_len(n) - 0.44) / (n + 0.12)
  obs <- data.frame(T = 1 / pp, q = sort(x))
  p <- ggplot2::ggplot()
  if (nrow(ci)) {
    band <- data.frame(T = ci$T, lower = ci$lower, upper = ci$upper)
    p <- p + ggplot2::geom_ribbon(data = band, ggplot2::aes(T, ymin = lower, ymax = upper),
                                  fill = PALETTE[1], alpha = 0.15)
  }
  p +
    ggplot2::geom_hline(yintercept = 0, colour = PALETTE[3], linetype = "dotted") +
    ggplot2::geom_line(data = curves, ggplot2::aes(T, q, colour = distribution,
                                                   linewidth = distribution == selected)) +
    ggplot2::geom_point(data = obs, ggplot2::aes(T, q), size = 1.3) +
    ggplot2::scale_x_log10(breaks = c(2, 5, 10, 20, 50, 100)) +
    ggplot2::scale_colour_manual(values = PALETTE) +
    ggplot2::scale_linewidth_manual(values = c(0.6, 1.3), guide = "none") +
    ggplot2::labs(title = title, x = "Return period T (years)",
                  y = sprintf("T-year low flow (%s)", unit),
                  subtitle = sprintf("band: 95%% interval on %s (selected by AIC); lower = more severe drought",
                                     selected)) +
    theme_lowflow()
}

#' Capture rate per candidate window, with the chosen one highlighted.
fig_season_scan <- function(scan, chosen = NULL, title = "") {
  s <- utils::head(scan[order(-scan$capture_rate, scan$width), ], 12)
  s$window <- factor(s$window, levels = rev(s$window))
  is_chosen <- if (is.null(chosen)) rep(FALSE, nrow(s)) else as.character(s$window) == chosen
  s$flag <- ifelse(is_chosen, "chosen",
                   ifelse(s$capture_rate >= 0.8, "adequate", "too low"))
  ggplot2::ggplot(s, ggplot2::aes(capture_rate, window, fill = flag)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_vline(xintercept = 0.8, linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.0f%%", 100 * capture_rate)),
                       hjust = -0.15, size = 3.2) +
    ggplot2::scale_fill_manual(values = c(chosen = PALETTE[3], adequate = PALETTE[2],
                                          `too low` = PALETTE[6])) +
    ggplot2::scale_x_continuous(labels = function(v) sprintf("%.0f%%", 100 * v),
                                limits = c(0, 1.12)) +
    ggplot2::labs(title = title, x = "Annual minima falling inside the window", y = NULL,
                  subtitle = "dashed line: the 80% floor below which the window is not describing the drought") +
    theme_lowflow() + ggplot2::theme(legend.position = "bottom")
}

#' Droughtiness pseudo-observations, with the model-free dependence measures.
fig_pseudo_obs <- function(u, name_x, name_y, tau, lambda, title = "") {
  d <- data.frame(u = u[, 1], v = u[, 2])
  ggplot2::ggplot(d, ggplot2::aes(u, v)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey50") +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::annotate("text", x = 0.02, y = 0.97, hjust = 0, size = 3.4,
                      label = sprintf("tau = %.3f\nlambda_drought (CFG) = %.3f", tau, lambda)) +
    ggplot2::labs(title = title, x = sprintf("%s droughtiness u", name_x),
                  y = sprintf("%s droughtiness v", name_y),
                  subtitle = "upper right = both rivers in drought") +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_lowflow()
}

#' Joint AND return-level contours in discharge space.
fig_and_curves <- function(curves, design, mld, obs, name_x, name_y,
                           title = "", unit = "m^3/s", copula_name = "") {
  curves$Tf <- factor(curves$T)
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(data = obs, ggplot2::aes(x, y), size = 1.5, colour = "black") +
    ggplot2::geom_path(data = curves, ggplot2::aes(x, y, colour = Tf, group = Tf),
                       linewidth = 0.7) +
    ggplot2::geom_point(data = design, ggplot2::aes(x, y), shape = 22, fill = NA,
                        colour = PALETTE[3], size = 2.4, stroke = 0.9)
  if (!is.null(mld) && nrow(mld)) {
    p <- p + ggplot2::geom_point(data = mld, ggplot2::aes(x, y), shape = 23,
                                 fill = PALETTE[3], colour = PALETTE[3], size = 2.2)
  }
  p + ggplot2::scale_colour_viridis_d(name = "T (yr)", end = 0.92) +
    ggplot2::labs(
      title = sprintf("%s%s", title, if (nzchar(copula_name))
        sprintf(": joint AND return levels (%s copula)", copula_name) else ""),
      subtitle = sprintf("P(Q_%s < q_%s and Q_%s < q_%s) = 1/T   -   toward the origin = both rivers drier",
                         name_x, name_x, name_y, name_y),
      x = sprintf("%s discharge (%s)", name_x, unit),
      y = sprintf("%s discharge (%s)", name_y, unit),
      caption = "square = equal-droughtiness point (u = v); diamond = most likely realisation on the contour") +
    theme_lowflow()
}

#' Save a figure as PNG and PDF, and record both paths.
save_fig <- function(p, outdir, stem, width = 7.5, height = 4.6) {
  png <- file.path(outdir, paste0(stem, ".png"))
  pdf <- file.path(outdir, paste0(stem, ".pdf"))
  suppressWarnings({
    ggplot2::ggsave(png, plot = p, width = width, height = height, dpi = 200)
    ggplot2::ggsave(pdf, plot = p, width = width, height = height)
  })
  c(png, pdf)
}
