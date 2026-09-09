## ---------------------------------------------------------------------------
## 07_exhibits.R
## Assembles the LaTeX table bodies used in the manuscript from the outputs of
## steps 04 and 05, so that no number is typed by hand into the .tex file.
## Inputs:  data/derived/boundary_estimates.csv, output/tables/*.csv
## Outputs: output/tables/table1_body.tex, table3_body.tex, tableC1_body.tex
## ---------------------------------------------------------------------------

est <- read_csv(file.path(DIR_DERIVED, "boundary_estimates.csv"), show_col_types = FALSE)
cnt <- file.path(DIR_DERIVED, "boundary_estimates_with_counts.csv")
if (file.exists(cnt)) {
  est <- left_join(est, read_csv(cnt, show_col_types = FALSE) |>
                     select(boundary_id, n_treated, n_control), by = "boundary_id")
}
est <- est |> mutate(n_eff = n_eff_left + n_eff_right) |> arrange(tau)

stars <- function(p) ifelse(p < 0.01, "\\sym{***}",
                     ifelse(p < 0.05, "\\sym{**}",
                     ifelse(p < 0.10, "\\sym{*}", "")))
fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)

## --- Table 1, Panel A ---------------------------------------------------------
body1 <- sprintf("%s & $%s$%s & (%s) & %s & [$%s$, $%s$] & %s & %s \\\\",
                 gsub("_", "--", est$boundary_id),
                 fmt(est$tau), stars(est$pval), fmt(est$se), fmt(est$pval),
                 fmt(est$ci_lower), fmt(est$ci_upper),
                 formatC(est$n_crimes, big.mark = ",", format = "d"),
                 formatC(est$n_eff,    big.mark = ",", format = "d"))
writeLines(body1, file.path(DIR_OUT_TAB, "table1_body.tex"))

## --- Appendix Table C.1 -------------------------------------------------------
summ <- function(x) c(mean(x), sd(x), min(x), median(x), max(x))
rows <- list(
  "Incidents within 2,000 m of the boundary"   = summ(est$n_crimes),
  "Incidents within the MSE-optimal bandwidth" = summ(est$n_treated + est$n_control),
  "\\quad of which treated"                    = summ(est$n_treated),
  "\\quad of which control"                    = summ(est$n_control),
  "Share treated within bandwidth (\\%)"       = summ(100*est$n_treated/(est$n_treated+est$n_control)),
  "Effective observations, bias-corrected"     = summ(est$n_eff),
  "MSE-optimal bandwidth (m)"                  = summ(est$bandwidth),
  "Point estimate ($\\hat{\\tau}_b$)"          = summ(est$tau),
  "Standard error"                             = summ(est$se),
  "Absolute effect size"                       = summ(abs(est$tau))
)
bodyC <- vapply(names(rows), function(nm) {
  v <- rows[[nm]]
  d <- if (max(abs(v)) > 20) 0 else 3
  sprintf("%s & %s & %s & %s & %s & %s & 15 \\\\", nm,
          fmt(v[1], d), fmt(v[2], d), fmt(v[3], d), fmt(v[4], d), fmt(v[5], d))
}, character(1))
writeLines(unname(bodyC), file.path(DIR_OUT_TAB, "tableC1_body.tex"))

cat("\nLaTeX table bodies written to output/tables/.\n")
cat("Totals for cross-checking against the manuscript:\n")
cat("  incidents within 2 km        :", sum(est$n_crimes), " (paper: 20,729)\n")
if (!is.null(est$n_treated))
  cat("  incidents within bandwidth   :", sum(est$n_treated + est$n_control), " (paper: 14,260)\n")
cat("  effective observations       :", sum(est$n_eff), " (paper: 4,451)\n")
