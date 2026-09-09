## ---------------------------------------------------------------------------
## _config.R -- paths, constants and package loading
## Replication package for:
##   Galvis Arias, N. "Does Urban Water Scarcity Increase Domestic Violence?
##   Evidence from Water Rationing in Bogota."
##
## EDIT ONLY THE TWO LINES MARKED "USER SETTING".
## ---------------------------------------------------------------------------

## USER SETTING 1: absolute path to the root of this replication package.
PROJECT_ROOT <- normalizePath(".", mustWork = TRUE)

## USER SETTING 2: absolute path to the restricted incident-level extract.
## Leave as NA to run only the steps that do not require the restricted data
## (steps 01, 04, 05 and 07 reproduce every reported table and figure from the
## derived files shipped in data/derived/).  See docs/data_access.md.
RESTRICTED_INCIDENTS <- NA_character_

## --- Derived paths (do not edit) --------------------------------------------
DIR_PUBLIC    <- file.path(PROJECT_ROOT, "data", "public")
DIR_DERIVED   <- file.path(PROJECT_ROOT, "data", "derived")
DIR_OUT_TAB   <- file.path(PROJECT_ROOT, "output", "tables")
DIR_OUT_FIG   <- file.path(PROJECT_ROOT, "output", "figures")
for (d in c(DIR_OUT_TAB, DIR_OUT_FIG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

## --- Analysis constants ------------------------------------------------------
## These are the exact values used in the paper. Changing them changes the
## reported numbers; they are collected here so that no magic number appears
## inside an estimation script.
CRS_METRIC       <- 3116     # EPSG:3116, MAGNA-SIRGAS / Colombia Bogota zone (metres)
BIN_WIDTH_M      <- 50       # signed-distance bin width, Section 4.5.1
MAX_DIST_M       <- 2000     # distance screen, Section 3.4 step 2
ADJACENCY_TOL_M  <- 1        # GIS tolerance for shared-boundary detection, Appendix A
DONUT_RADIUS_M   <- 100      # donut RD exclusion window, Table 3 column (3)
MIN_BOUNDARY_N   <- 500      # "large boundary" screen, Table 3 column (2)
BASELINE_RATE    <- 2.03     # control-side incidents per boundary-day, Section 5.2
POWER_MULTIPLIER <- 2.80     # (z_{0.975} + z_{0.80}) for the MDE, Section 5.2
STUDY_START      <- as.Date("2024-04-11")
STUDY_END        <- as.Date("2025-04-01")

## --- Packages ----------------------------------------------------------------
## Required on BOTH execution paths.
.pkgs_core <- c("dplyr", "tidyr", "readr", "purrr", "tibble", "ggplot2", "sf")
## Required only on the FULL path, which re-estimates from the restricted
## microdata. The derived-data path reproduces every reported table and figure
## from the shipped boundary estimates and never calls these, so it must not be
## blocked by their absence.
.pkgs_full <- c("rdrobust", "rddensity")

.missing_core <- .pkgs_core[!vapply(.pkgs_core, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing_core)) {
  stop("Missing R packages required on every path: ", paste(.missing_core, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste0('"', .missing_core, '"', collapse = ", "), "))")
}
.missing_full <- .pkgs_full[!vapply(.pkgs_full, requireNamespace, logical(1), quietly = TRUE)]
if (length(.missing_full) && !is.na(RESTRICTED_INCIDENTS)) {
  stop("The full path needs: ", paste(.missing_full, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste0('"', .missing_full, '"', collapse = ", "), "))",
       "\nOr leave RESTRICTED_INCIDENTS = NA to run the derived-data path.")
}
if (length(.missing_full))
  message("Note: ", paste(.missing_full, collapse = ", "),
          " not installed. Derived-data path only; steps 02-03 unavailable.")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

message("Config loaded. Restricted data: ",
        if (is.na(RESTRICTED_INCIDENTS)) "NOT SUPPLIED (derived-data path)" else RESTRICTED_INCIDENTS)
