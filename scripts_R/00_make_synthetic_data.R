# Generate three synthetic daily-discharge workbooks in data/.
#
# These stand in for your own records so the pipeline runs end to end out of the
# box: two gauges in a shared basin (a large river and a smaller tributary) plus
# a transient projection series for script 03. The series carry the features the
# pipeline is built to survive: a seasonal cycle with late-summer minima, AR(1)
# persistence, a climate signal shared between the two gauges, a multi-month
# outage, a year with only three months of record, scattered missing days and a
# handful of zero-flow days.
#
# Replace data/*.xlsx with your own records -- same two columns, date and flow
# -- and point the scripts at them. Nothing downstream is specific to these.
#
# Run from the repository root:  Rscript scripts_R/00_make_synthetic_data.R

suppressPackageStartupMessages(library(writexl))

SEED <- 20240617
set.seed(SEED)

seasonal_log <- function(doy, lo_doy, hi, lo) {
  ph <- 2 * pi * (doy - lo_doy) / 365.25
  log(lo) + (log(hi) - log(lo)) * 0.5 * (1 - cos(ph))
}

ar1 <- function(n, rho, sd) {
  e <- stats::rnorm(n, 0, sd * sqrt(1 - rho^2))
  x <- numeric(n)
  x[1] <- stats::rnorm(1, 0, sd)
  for (i in 2:n) x[i] <- rho * x[i - 1] + e[i]
  x
}

make_river <- function(dates, hi, lo, lo_doy, sd, rho, shared, share_w) {
  doy <- as.integer(format(dates, "%j"))
  base <- seasonal_log(doy, lo_doy, hi, lo)
  idio <- ar1(length(dates), rho, sd)
  ann <- shared[as.character(as.integer(format(dates, "%Y")))]
  exp(base + share_w * ann + sqrt(max(0, 1 - share_w^2)) * idio)
}

outdir <- "data"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

dates <- seq(as.Date("1970-01-01"), as.Date("2024-12-31"), by = "day")
years <- sort(unique(as.integer(format(dates, "%Y"))))
shared <- stats::setNames(ar1(length(years), 0.25, 0.38), as.character(years))

big <- make_river(dates, 2600, 480, 213, 0.34, 0.93, shared, 0.75)
small <- make_river(dates, 95, 4.2, 220, 0.52, 0.90, shared, 0.62)

a <- data.frame(date = dates, flow = round(big, 1))
b <- data.frame(date = dates, flow = round(small, 2))

# a multi-month outage at the smaller gauge, and a year with only three months
b$flow[dates >= as.Date("1983-04-10") & dates <= as.Date("1983-11-02")] <- NA
a$flow[as.integer(format(dates, "%Y")) == 1991 &
         as.integer(format(dates, "%m")) > 3] <- NA
# scattered missing days
idx <- which(!is.na(b$flow))
b$flow[sample(idx, 140)] <- NA
# a handful of recorded zeros, which the pipeline keeps by default
dry <- which(!is.na(b$flow) & b$flow < 1)
if (length(dry)) b$flow[utils::head(dry, 6)] <- 0

pdates <- seq(as.Date("2015-01-01"), as.Date("2100-12-31"), by = "day")
pyears <- sort(unique(as.integer(format(pdates, "%Y"))))
trend <- stats::setNames(seq(0, -0.55, length.out = length(pyears)) +
                           ar1(length(pyears), 0.2, 0.33), as.character(pyears))
proj <- make_river(pdates, 90, 3.9, 224, 0.55, 0.90, trend, 0.70)
p <- data.frame(date = pdates, flow = round(proj, 2))

for (item in list(list(a, "station_a"), list(b, "station_b"), list(p, "projection"))) {
  df <- item[[1]]; nm <- item[[2]]
  writexl::write_xlsx(df, file.path(outdir, paste0(nm, ".xlsx")))
  cat(sprintf("wrote %s/%s.xlsx  (%d rows, %d missing)\n",
              outdir, nm, nrow(df), sum(is.na(df$flow))))
}
