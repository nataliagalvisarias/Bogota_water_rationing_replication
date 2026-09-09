## ---------------------------------------------------------------------------
## 00_master.R -- run the full replication
##
## Usage, from the package root:
##     Rscript code/00_master.R
##
## Two execution paths:
##   (A) DERIVED-DATA PATH (default).  RESTRICTED_INCIDENTS is NA in _config.R.
##       Runs steps 01, 04, 05, 07.  Reproduces every table and figure reported
##       in the paper from the derived files in data/derived/.  Requires no
##       restricted data.  Expected runtime: under 2 minutes.
##   (B) FULL PATH.  Set RESTRICTED_INCIDENTS in _config.R to the incident-level
##       extract obtained from the Bogota Secretariat of Security (see
##       docs/data_access.md).  Runs steps 01-07, rebuilding the derived files
##       from the raw records.  Expected runtime: 20-40 minutes.
##
## Every step writes its outputs to output/ and prints the reported value
## alongside the recomputed value, so that any divergence is visible
## immediately rather than discovered later.
## ---------------------------------------------------------------------------

source("code/_config.R")

full_path <- !is.na(RESTRICTED_INCIDENTS)

run <- function(f) { message("\n>>> ", f); source(file.path("code", f), echo = FALSE) }

run("01_build_adjacency.R")            # public data only
if (full_path) {
  run("02_build_panel.R")              # requires restricted data
  run("03_estimate_boundary_rd.R")     # requires the panel from step 02
} else {
  message("\n[skipped] 02_build_panel.R and 03_estimate_boundary_rd.R ",
          "(restricted data not supplied); using shipped derived estimates.")
}
run("04_aggregate_inference.R")        # Table 1 Panels B-C, Table 4, Figs 3-5
run("05_robustness.R")                 # Table 3
run("07_exhibits.R")                   # assembles output/tables and output/figures

## Final submission estimates: the balanced, zero-filled, deduplicated panel that
## produces Tables 1, 2, 3 and C.1 and Figures 3-7 in the manuscript. These steps
## need no restricted data and run on both paths.
message("\n>>> 10_final_submission_estimates.R")
local({ Sys.setenv(ROOT = PROJECT_ROOT); source("code/10_final_submission_estimates.R") })
message("\n>>> 11_final_figures_and_checks.R")
local({ Sys.setenv(ROOT = PROJECT_ROOT); source("code/11_final_figures_and_checks.R") })

## Exact design-based inference. This is the PRIMARY inferential result reported
## in Sections 4.7 and 5.5 and in Table 4; the cluster-robust and date-block
## bootstrap standard errors are retained only as points of comparison.
message("\n>>> 12_randomization_inference.R")
local({ Sys.setenv(ROOT = PROJECT_ROOT); source("code/12_randomization_inference.R") })

message("\nThe legacy Python design checks in 06_design_checks.py are superseded by")
message("step 12 and are retained only for the balance and density diagnostics.")
message("\nDone. Outputs in output/tables and output/figures.")
