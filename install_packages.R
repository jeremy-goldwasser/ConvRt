## ============================================================================
## install_packages.R
##
## Install everything needed to use the ConvRt package and run examples/demo.R.
## ConvRt is an R package, so the simplest route is to let install_github pull
## the dependencies declared in DESCRIPTION automatically:
##
##   Rscript -e 'remotes::install_github("jeremy-goldwasser/ConvRt")'
##
## This script does the same thing explicitly (handy when developing locally),
## and additionally installs the few packages the demo uses for plotting.
##
## NOTE: this is a software-package repo, not a benchmark study -- so it does NOT
## install the other Rt estimators (EpiEstim, rtestim, EpiLPS, EpiNow2, rstan,
## estimateR). None of those are needed to fit or demo ConvRt.
## ============================================================================

cat("Installing R package dependencies for ConvRt...\n\n")

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s from CRAN...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat(sprintf("  %s already installed (v%s)\n", pkg, packageVersion(pkg)))
  }
}

# ---- ConvRt runtime dependencies (DESCRIPTION: Depends + Imports) ----
#   splines, stats, utils are base R (shipped with R) -- nothing to install.
cran_pkgs <- c(
  "dplyr",      # data manipulation        (Depends)
  "tibble",     # data frames              (Depends)
  "tidyr",      # data manipulation        (Depends)
  # ---- demo / optional features ----
  "ggplot2",    # plotting (examples/demo.R)
  "patchwork",  # stacking demo panels    (examples/demo.R)
  "mgcv",       # OPTIONAL: gi_deconvolve_exposures(link = "log")
  "remotes"     # to install ConvRt itself from GitHub
)
# NOTE: the trend-filter variant (fit_convrt_tf_*) needs CVXR, which is NOT
# installed here -- it is a heavy, rarely-used optional dependency.  Those
# functions print an install hint (install.packages("CVXR")) if you call them
# without it; everything else in ConvRt works without CVXR.

cat("--- CRAN packages ---\n")
for (pkg in cran_pkgs) install_if_missing(pkg)

# ---- ConvRt itself ----
cat("\n--- ConvRt ---\n")
if (!requireNamespace("ConvRt", quietly = TRUE)) {
  cat("  Installing ConvRt from GitHub...\n")
  remotes::install_github("jeremy-goldwasser/ConvRt")
} else {
  cat(sprintf("  ConvRt already installed (v%s)\n", packageVersion("ConvRt")))
}

# ---- Verify ----
cat("\n--- Verification ---\n")
ok <- TRUE
for (pkg in c("dplyr", "tibble", "tidyr", "splines", "ggplot2", "patchwork", "ConvRt")) {
  loaded <- tryCatch({
    suppressPackageStartupMessages(library(pkg, character.only = TRUE)); TRUE
  }, error = function(e) FALSE)
  cat(sprintf("  %-12s  %s\n", pkg, if (loaded) "OK" else "FAILED"))
  if (!loaded) ok <- FALSE
}
cat(if (ok) "\nAll set.\n" else "\nWARNING: some packages failed to load.\n")
