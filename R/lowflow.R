# Loader. Source this one file and the whole interface is available:
#
#   source("R/lowflow.R")
#
# Run from the repository root so that relative data/ and out/ paths resolve.
#
# Needs: readxl, zoo, lmomco, fitdistrplus, copula, Kendall, ggplot2
# (writexl as well, for scripts_R/00_make_synthetic_data.R). Install with:
#
#   install.packages(c("readxl", "writexl", "zoo", "lmomco", "fitdistrplus",
#                      "copula", "Kendall", "ggplot2"))

LOWFLOW_VERSION <- "1.2.0"

.lowflow_dir <- local({
  d <- NULL
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) { d <- dirname(normalizePath(of, mustWork = FALSE)); break }
  }
  if (is.null(d) || !file.exists(file.path(d, "dataio.R"))) {
    cand <- file.path(getwd(), "R")
    if (file.exists(file.path(cand, "dataio.R"))) d <- cand
  }
  if (is.null(d) || !file.exists(file.path(d, "dataio.R"))) {
    stop("could not locate the R/ directory; run from the repository root, ",
         "or source R/lowflow.R by its full path")
  }
  d
})

for (.f in c("dataio.R", "marginals.R", "copulas.R", "figures.R", "pipeline.R")) {
  source(file.path(.lowflow_dir, .f))
}
rm(.f)

#' Report the loaded version and whether every dependency is present.
lowflow_check <- function() {
  need <- c("readxl", "zoo", "lmomco", "fitdistrplus", "copula",
            "Kendall", "ggplot2")
  have <- vapply(need, requireNamespace, logical(1), quietly = TRUE)
  cat(sprintf("lowflow (R) %s loaded from %s\n", LOWFLOW_VERSION, .lowflow_dir))
  for (i in seq_along(need)) {
    cat(sprintf("  %-14s %s\n", need[i], if (have[i]) "ok" else "MISSING"))
  }
  if (!all(have)) {
    cat("\ninstall.packages(c(", paste(sprintf('"%s"', need[!have]), collapse = ", "), "))\n")
  }
  invisible(all(have))
}
