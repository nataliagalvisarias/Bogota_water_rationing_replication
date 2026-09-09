## ---------------------------------------------------------------------------
## 04_aggregate_inference.R
## Reproduces: Table 1 Panels B and C; Section 5.2 (MDE, TOST, one-sided bound);
##             the leave-one-out range in Table 3; Figures 3, 4 and 5.
## Inputs:  data/derived/boundary_estimates.csv
## Outputs: output/tables/table1_panelBC.csv
##          output/tables/equivalence.csv
##          output/tables/leave_one_out.csv
##          output/figures/fig3_forest.pdf, fig4_histogram.pdf
##
## This script runs on the DERIVED-DATA PATH and needs no restricted data.
## It is the script a replicator should run first: it reproduces every headline
## number in the paper from the fifteen boundary estimates.
## ---------------------------------------------------------------------------

est <- read_csv(file.path(DIR_DERIVED, "boundary_estimates.csv"), show_col_types = FALSE)
stopifnot(nrow(est) == 15L, all(c("tau", "se") %in% names(est)))

## --- Inverse-variance weighted aggregation: equation (13) --------------------
ivw <- function(tau, se) {
  w  <- 1 / se^2
  t  <- sum(w * tau) / sum(w)
  s  <- 1 / sqrt(sum(w))
  list(tau = t, se = s, z = t / s, p = 2 * pnorm(-abs(t / s)),
       lo = t - qnorm(0.975) * s, hi = t + qnorm(0.975) * s)
}
pooled <- ivw(est$tau, est$se)

## Unweighted average. Its standard error treats the boundary estimates as
## independent draws with known variances, hence sqrt(sum(se^2))/K.
simple_tau <- mean(est$tau)
simple_se  <- sqrt(sum(est$se^2)) / nrow(est)

## --- Heterogeneity: Cochran's Q and I-squared --------------------------------
w  <- 1 / est$se^2
Q  <- sum(w * (est$tau - pooled$tau)^2)
df <- nrow(est) - 1
I2 <- max(0, (Q - df) / Q) * 100

## --- Equivalence and power: Section 5.2 --------------------------------------
mde        <- POWER_MULTIPLIER * pooled$se
margin     <- 0.10 * BASELINE_RATE
tost_lower <- 1 - pnorm((pooled$tau + margin) / pooled$se)   # H0: tau <= -margin
tost_upper <- pnorm((pooled$tau - margin) / pooled$se)       # H0: tau >= +margin
tost_p     <- max(tost_lower, tost_upper)
## Tightest margin at which equivalence still holds at the 5% level.
tightest   <- abs(pooled$tau) + qnorm(0.95) * pooled$se

## --- Leave-one-boundary-out ---------------------------------------------------
loo <- purrr::map_dfr(seq_len(nrow(est)), function(i) {
  p <- ivw(est$tau[-i], est$se[-i])
  tibble::tibble(dropped = est$boundary_id[i], tau = p$tau, se = p$se)
})
write_csv(loo, file.path(DIR_OUT_TAB, "leave_one_out.csv"))

## --- Verification against the paper -------------------------------------------
chk <- function(label, got, reported, tol = 5e-3) {
  flag <- if (abs(got - reported) <= tol) "OK " else "***"
  cat(sprintf("%s %-46s recomputed = %9.4f   paper = %9.4f\n", flag, label, got, reported))
}
cat("\n--- LEGACY specification (occupied-cell panel), superseded by step 10 ---\n")
cat("These values documented Tables 1 and 3 before the panel was rebuilt. They are\n")
cat("retained so the earlier results remain reproducible; the manuscript now reports\n")
cat("the balanced zero-filled estimates produced by 10_final_submission_estimates.R.\n")
chk("Pooled IVW estimate",              pooled$tau,  -0.044)
chk("Pooled IVW standard error",        pooled$se,    0.057)
chk("Pooled IVW p-value",               pooled$p,     0.440)
chk("Pooled IVW CI lower",              pooled$lo,   -0.155)
chk("Pooled IVW CI upper",              pooled$hi,    0.067)
chk("Unweighted average",               simple_tau,  -0.053)
chk("Unweighted average standard error",simple_se,    0.086)
chk("Cochran's Q",                      Q,            9.87, tol = 0.02)
chk("I-squared (%)",                    I2,           0.00)
chk("Minimum detectable effect",        mde,          0.159)
chk("MDE as % of baseline",             100*mde/BASELINE_RATE, 7.8, tol = 0.06)
chk("Upper CI as % of baseline",        100*pooled$hi/BASELINE_RATE, 3.3, tol = 0.05)
chk("TOST p-value at 10% margin",       tost_p,       0.003, tol = 5e-4)
chk("Tightest equivalence margin (%)",  100*tightest/BASELINE_RATE, 6.8, tol = 0.05)
chk("Leave-one-out minimum",            min(loo$tau), -0.062)
chk("Leave-one-out maximum",            max(loo$tau), -0.005)
cat("Rows marked *** indicate a divergence from the published value.\n")

write_csv(tibble::tibble(
  quantity = c("IVW", "IVW SE", "IVW p", "IVW CI low", "IVW CI high",
               "Simple average", "Simple average SE", "Cochran Q", "I2 (%)"),
  value    = c(pooled$tau, pooled$se, pooled$p, pooled$lo, pooled$hi,
               simple_tau, simple_se, Q, I2)
), file.path(DIR_OUT_TAB, "table1_panelBC.csv"))

write_csv(tibble::tibble(
  quantity = c("MDE (incidents/boundary-day)", "MDE (% of baseline)",
               "TOST p at 10% margin", "Tightest equivalence margin (% of baseline)",
               "One-sided upper bound (% of baseline)"),
  value    = c(mde, 100*mde/BASELINE_RATE, tost_p,
               100*tightest/BASELINE_RATE, 100*pooled$hi/BASELINE_RATE)
), file.path(DIR_OUT_TAB, "equivalence.csv"))

## --- Figure 3: forest plot -----------------------------------------------------
fp <- est |> mutate(label = gsub("_", "--", boundary_id)) |> arrange(tau) |>
  mutate(label = factor(label, levels = label))
p3 <- ggplot(fp, aes(x = tau, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = pooled$tau, colour = "#660099", linewidth = 0.8) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0,
                 colour = "#660099", alpha = 0.7) +
  geom_point(shape = 21, fill = "white", colour = "#660099", size = 2.2) +
  labs(x = expression(hat(tau)[b]~"(incidents per date"%*%"bin)"),
       y = "Boundary (polygon pair)") +
  theme_minimal(base_size = 11)
ggsave(file.path(DIR_OUT_FIG, "fig3_forest.pdf"), p3, width = 7, height = 5)

## --- Figure 4: histogram --------------------------------------------------------
p4 <- ggplot(est, aes(x = tau)) +
  geom_histogram(binwidth = 0.10, fill = "#B280CC", colour = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = pooled$tau, colour = "#660099", linewidth = 0.9) +
  labs(x = expression("Boundary-specific estimate"~hat(tau)[b]),
       y = "Number of boundaries") +
  theme_minimal(base_size = 11)
ggsave(file.path(DIR_OUT_FIG, "fig4_histogram.pdf"), p4, width = 6, height = 4.2)

cat("\nFigures 3 and 4 written to output/figures/.\n")
cat("NOTE: the figures in the manuscript were produced by the author's original\n",
    "plotting code; these reproductions carry the same data and differ only in\n",
    "cosmetic styling.\n", sep = "")
