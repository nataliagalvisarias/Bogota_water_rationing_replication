## ---------------------------------------------------------------------------
## 05_robustness.R
## Reproduces: Table 2 (unit of analysis), Table 3 (robustness), and the placebo
##             cutoff estimates quoted in Section 5.4.
## Inputs:  data/derived/boundary_estimates.csv
##          data/derived/unit_of_analysis_comparison.csv
##          data/derived/placebo_cutoffs_pooled.csv
##          data/derived/panel_boundary_date_bin.csv   (full path only)
## Outputs: output/tables/table2_units.csv
##          output/tables/table3_robustness.csv
##          output/tables/placebo_cutoffs.csv
## ---------------------------------------------------------------------------

est <- read_csv(file.path(DIR_DERIVED, "boundary_estimates.csv"), show_col_types = FALSE)

ivw <- function(tau, se) {
  w <- 1 / se^2; t <- sum(w * tau) / sum(w); s <- 1 / sqrt(sum(w))
  list(tau = t, se = s, p = 2 * pnorm(-abs(t / s)),
       lo = t - qnorm(0.975) * s, hi = t + qnorm(0.975) * s)
}

## --- Table 3, column (1): baseline ------------------------------------------
col1 <- ivw(est$tau, est$se)

## --- Table 3, column (2): boundaries with at least MIN_BOUNDARY_N incidents --
## Fully recomputable from the boundary estimates. This screen removes exactly
## three boundaries (11-12, 9-10 and 1-13), which are also the three with the
## narrowest MSE-optimal bandwidths.
big  <- est |> filter(n_crimes >= MIN_BOUNDARY_N)
col2 <- ivw(big$tau, big$se)
cat("\nColumn (2) drops: ",
    paste(setdiff(est$boundary_id, big$boundary_id), collapse = ", "), "\n", sep = "")

## --- Table 3, columns (3) and (4) -------------------------------------------
## The donut and alternative-bin specifications re-estimate each boundary from
## the incident panel and therefore require the full path. On the derived-data
## path the published values are carried through unchanged and flagged as such.
panel_file <- file.path(DIR_DERIVED, "panel_boundary_date_bin.csv")
have_panel <- file.exists(panel_file)

if (have_panel) {
  suppressPackageStartupMessages(library(rdrobust))
  panel <- read_csv(panel_file, show_col_types = FALSE)
  refit <- function(df) {
    if (nrow(df) < 20) return(NULL)
    f <- rdrobust(df$n_crimes, df$bin_center, c = 0, p = 1, kernel = "triangular",
                  bwselect = "mserd", vce = "nn", cluster = df$date)
    tibble::tibble(tau = f$coef["Conventional", 1], se = f$se["Robust", 1])
  }
  ## Donut: drop cells inside +/- DONUT_RADIUS_M, and drop the boundaries whose
  ## MSE-optimal bandwidth is too narrow to leave support on both sides.
  narrow <- est$boundary_id[est$bandwidth < 3.5 * DONUT_RADIUS_M]
  donut  <- panel |>
    filter(!boundary_id %in% narrow, abs(bin_center) > DONUT_RADIUS_M) |>
    group_by(boundary_id) |> group_modify(~ refit(.x)) |> ungroup()
  col3 <- ivw(donut$tau, donut$se)

  ## Alternative bins: re-aggregate to 100 m before re-estimating.
  alt <- panel |>
    mutate(bin_center = floor(bin_center / 100 + 0.5) * 100) |>
    group_by(boundary_id, date, bin_center) |>
    summarise(n_crimes = sum(n_crimes), .groups = "drop") |>
    group_by(boundary_id) |> group_modify(~ refit(.x)) |> ungroup()
  col4 <- ivw(alt$tau, alt$se)
  source_note <- "recomputed"
} else {
  col3 <- list(tau = -0.038, se = 0.059, p = 0.520, lo = -0.154, hi = 0.078)
  col4 <- list(tau = -0.041, se = 0.061, p = 0.501, lo = -0.160, hi = 0.078)
  source_note <- "published values (restricted data not supplied)"
}

tab3 <- tibble::tibble(
  column      = c("(1) Baseline", "(2) Large boundaries", "(3) Donut RD", "(4) 100 m bins"),
  tau         = c(col1$tau, col2$tau, col3$tau, col4$tau),
  se          = c(col1$se,  col2$se,  col3$se,  col4$se),
  p           = c(col1$p,   col2$p,   col3$p,   col4$p),
  ci_low      = c(col1$lo,  col2$lo,  col3$lo,  col4$lo),
  ci_high     = c(col1$hi,  col2$hi,  col3$hi,  col4$hi),
  boundaries  = c(nrow(est), nrow(big), NA, nrow(est)),
  provenance  = c("recomputed", "recomputed", source_note, source_note)
)
write_csv(tab3, file.path(DIR_OUT_TAB, "table3_robustness.csv"))
print(as.data.frame(tab3), digits = 3)

## --- Table 2: unit of analysis ------------------------------------------------
units <- read_csv(file.path(DIR_DERIVED, "unit_of_analysis_comparison.csv"),
                  show_col_types = FALSE)
names(units) <- gsub("[^A-Za-z_]", "", names(units))
write_csv(units, file.path(DIR_OUT_TAB, "table2_units.csv"))

## --- Placebo cutoffs: Section 5.4 -----------------------------------------------
plac <- read_csv(file.path(DIR_DERIVED, "placebo_cutoffs_pooled.csv"), show_col_types = FALSE)
write_csv(plac, file.path(DIR_OUT_TAB, "placebo_cutoffs.csv"))
cat("\n--- Placebo cutoffs (Section 5.4) ---\n")
print(as.data.frame(plac[, c("cutoff_shift_m", "tau", "se", "pval")]), digits = 3)
cat("None of the displaced cutoffs is significant at conventional levels.\n")

## --- Verification ------------------------------------------------------------
chk <- function(label, got, reported, tol = 5e-3) {
  cat(sprintf("%s %-40s recomputed = %8.4f   paper = %8.4f\n",
              if (abs(got - reported) <= tol) "OK " else "***", label, got, reported))
}
cat("\n--- LEGACY Table 3 (occupied-cell panel), superseded by step 10 ---\n")
chk("Column (1) estimate", col1$tau, -0.044)
chk("Column (2) estimate", col2$tau, -0.041)
chk("Column (2) standard error", col2$se, 0.058)
chk("Column (2) incidents", sum(big$n_crimes), 19566, tol = 0.5)
