## ---------------------------------------------------------------------------
## 03_estimate_boundary_rd.R
## Reproduces: Table 1, Panel A (boundary-specific estimates) and the bandwidth
##             and effective-sample columns of Appendix Table C.1.
## Inputs:  data/derived/panel_boundary_date_bin.csv  (from step 02)
## Outputs: data/derived/boundary_estimates.csv       (overwrites the shipped file)
##
## RUNS ONLY ON THE FULL PATH.
## ---------------------------------------------------------------------------

if (is.na(RESTRICTED_INCIDENTS)) {
  stop("03_estimate_boundary_rd.R requires the panel built by 02_build_panel.R.")
}
suppressPackageStartupMessages(library(rdrobust))

panel <- read_csv(file.path(DIR_DERIVED, "panel_boundary_date_bin.csv"),
                  show_col_types = FALSE)

## Local linear RD with triangular kernel, MSE-optimal bandwidth and robust
## bias correction (Calonico, Cattaneo and Titiunik 2014). Standard errors use
## nearest-neighbour variance estimation clustered by date, which is what makes
## the within-date logic of the design operative in the inference
## (Section 4.5.1).
estimate_one <- function(df) {
  fit <- rdrobust(y       = df$n_crimes,
                  x       = df$bin_center,
                  c       = 0,
                  p       = 1,
                  kernel  = "triangular",
                  bwselect= "mserd",
                  vce     = "nn",
                  cluster = df$date)
  tibble::tibble(
    tau        = fit$coef["Conventional", 1],
    se         = fit$se["Robust", 1],
    pval       = fit$pv["Robust", 1],
    ci_lower   = fit$ci["Robust", 1],
    ci_upper   = fit$ci["Robust", 2],
    bandwidth  = fit$bws["h", 1],
    n_eff_left = fit$N_h[1],
    n_eff_right= fit$N_h[2]
  )
}

est <- panel |>
  group_by(boundary_id) |>
  group_modify(~ estimate_one(.x)) |>
  ungroup() |>
  left_join(count(panel, boundary_id, wt = n_crimes, name = "n_crimes"),
            by = "boundary_id") |>
  arrange(tau)

write_csv(est, file.path(DIR_DERIVED, "boundary_estimates.csv"))
print(as.data.frame(est), digits = 4)
cat("\n", nrow(est), " boundaries estimated (paper: 15).\n", sep = "")
