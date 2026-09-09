#!/usr/bin/env Rscript
## =============================================================================
## 10_final_submission_estimates.R
##
## Final estimation pipeline for
##   Galvis Arias, N. "Does Urban Water Scarcity Increase Domestic Violence?
##   Evidence from Water Rationing in Bogota."   (Journal of Urban Economics)
##
## Produces every number in Tables 1, 2, 3 and C.1 and every statistic quoted in
## Sections 3, 4 and 5, from one input file, in one pass, and writes the LaTeX
## table bodies directly so that no figure is ever typed by hand.
##
## -----------------------------------------------------------------------------
## THREE DESIGN DECISIONS, STATED UP FRONT
##
## 1. SIGN CONVENTION.  Equation (6) puts treated observations on the NEGATIVE
##    side of the running variable.  rdrobust's coefficient is the right limit
##    minus the left limit, i.e. CONTROL minus TREATED, which is the negative of
##    the treatment effect.  This script defines
##        tau = lim_{r->0-} E[Y|R=r] - lim_{r->0+} E[Y|R=r] = treated - control,
##    so that a POSITIVE tau means rationing raises reported violence.  Every
##    estimate below is on that convention.
##
## 2. STRICT, NON-OVERLAPPING BINS.  Bins are half-open and never span the
##    cutoff: the treated side is [-50,0), [-100,-50), ... and the control side
##    [0,50), [50,100), ..., with centres at -25, -75, ... and +25, +75, ...
##    Rounding to the nearest multiple of 50 would put a cell at the cutoff
##    holding observations from both sides.
##
## 3. ZERO-FILLED BALANCED PANEL.  The estimation panel is the complete
##    boundary x date x bin grid within +/-2,000 m, with zero-incidence cells
##    filled.  Counting only occupied cells conditions on a positive outcome and
##    attenuates the discontinuity toward zero.
##
## -----------------------------------------------------------------------------
## ESTIMATOR
##   Local linear with a triangular kernel at the MSE-optimal bandwidth, fitted
##   as one interacted weighted least squares per boundary,
##       y = a + tau*D + b1*x + b2*x*D,   D = 1{x < 0},   w = max(0, 1-|x|/h),
##   whose tau is algebraically identical to the difference of the two
##   side-specific intercepts.  Standard errors are cluster-robust by date,
##   which is the level at which the design's within-date comparison operates.
##
##   If the rdrobust package is installed it is used for MSE-optimal bandwidth
##   selection and its bandwidths are written to data/derived/bandwidths_mserd.csv.
##   If it is not, the shipped bandwidth file is read instead, so the script is
##   fully reproducible either way.  Set FORCE_PINNED_BW <- TRUE to use the
##   shipped bandwidths even when rdrobust is available.
## =============================================================================

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(readr)})
options(stringsAsFactors = FALSE)

ROOT   <- Sys.getenv("ROOT", ".")
DER    <- file.path(ROOT, "data", "derived")
TAB    <- file.path(ROOT, "tables")
OUT    <- file.path(ROOT, "output")
for (d in c(DER, TAB, OUT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

BIN_W        <- 50      # metres
MAX_D        <- 2000    # distance screen, metres
DONUT_R      <- 100     # donut exclusion radius, metres
MIN_BND_N    <- 500     # "large boundary" screen, incidents within 2,000 m
N_BOOT       <- 500
SEED         <- 20240411
EQUIV_MARGIN <- 0.10    # equivalence margin as a share of the control baseline
FORCE_PINNED_BW <- as.logical(Sys.getenv("FORCE_PINNED_BW", "FALSE"))

HAVE_RDROBUST <- requireNamespace("rdrobust", quietly = TRUE) && !FORCE_PINNED_BW
message("Bandwidth source: ", if (HAVE_RDROBUST) "rdrobust (MSE-optimal, recomputed)"
        else "data/derived/bandwidths_mserd.csv (pinned, produced by rdrobust)")

## -----------------------------------------------------------------------------
## 1. Panel construction
## -----------------------------------------------------------------------------
inc <- read_csv(file.path(DER, "rd_estimation_sample.csv"), show_col_types = FALSE) |>
  filter(in_estimation_sample == 1, is_duplicate_record == 0) |>
  mutate(date = as.Date(date))

N_INCIDENTS <- nrow(inc)
message(sprintf("Incidents in the deduplicated estimation sample: %s", format(N_INCIDENTS, big.mark = ",")))

## Strict half-open bins; no cell spans the cutoff.
bin_of <- function(x) sign(x) * (floor(abs(x) / BIN_W) * BIN_W + BIN_W / 2)
occ <- inc |> mutate(bin = bin_of(running_var)) |>
  count(boundary_id, date, bin, name = "n")

## Complete balanced grid within +/-MAX_D for every active boundary-date.
grid_bins <- c(seq(-MAX_D + BIN_W/2, -BIN_W/2, by = BIN_W),
               seq(BIN_W/2,  MAX_D - BIN_W/2, by = BIN_W))
panel <- occ |> distinct(boundary_id, date) |>
  tidyr::crossing(bin = grid_bins) |>
  left_join(occ, by = c("boundary_id", "date", "bin")) |>
  mutate(n = tidyr::replace_na(n, 0L))

N_CELLS <- nrow(panel)
PCT_ZERO <- 100 * mean(panel$n == 0)
BASELINE <- mean(panel$n[panel$bin > 0])          # control-side mean incidents per cell
message(sprintf("Balanced panel: %s cells across %s boundary-dates, %.1f%% zeros; control baseline %.4f incidents/cell",
                format(N_CELLS, big.mark = ","),
                format(n_distinct(paste(panel$boundary_id, panel$date)), big.mark = ","),
                PCT_ZERO, BASELINE))
write_csv(panel, file.path(DER, "panel_balanced_zerofilled.csv"))

## -----------------------------------------------------------------------------
## 2. Estimator
## -----------------------------------------------------------------------------

## Cluster-robust (CR1) sandwich variance for a weighted least squares fit.
cluster_vcov <- function(X, u, w, cl) {
  bread <- solve(crossprod(X * sqrt(w)))
  meat  <- matrix(0, ncol(X), ncol(X))
  for (g in split(seq_along(u), cl)) {
    s <- crossprod(X[g, , drop = FALSE], w[g] * u[g])
    meat <- meat + tcrossprod(s)
  }
  G <- length(unique(cl)); N <- length(u); K <- ncol(X)
  bread %*% meat %*% bread * (G / (G - 1)) * ((N - 1) / (N - K))
}

## Local linear RD at a fixed bandwidth.  tau is treated minus control.
rd_local_linear <- function(y, x, h, cl) {
  w <- pmax(0, 1 - abs(x) / h)
  k <- w > 0
  if (sum(k) < 8L) return(NULL)
  y <- y[k]; x <- x[k]; w <- w[k]; cl <- cl[k]
  D <- as.numeric(x < 0)                       # 1 on the treated side
  X <- cbind(`(Intercept)` = 1, treated = D, x = x, `x:treated` = x * D)
  if (qr(X * sqrt(w))$rank < ncol(X)) return(NULL)
  beta <- solve(crossprod(X * sqrt(w)), crossprod(X, w * y))
  u    <- as.vector(y - X %*% beta)
  V    <- cluster_vcov(X, u, w, cl)
  list(tau = unname(beta[2]), se = sqrt(V[2, 2]),
       n_eff = sum(k), n_cluster = length(unique(cl)),
       n_left = sum(x < 0), n_right = sum(x > 0),
       n_incidents_h = sum(y))
}

## MSE-optimal bandwidths: recomputed with rdrobust when available, else pinned.
bw_file <- file.path(DER, "bandwidths_mserd.csv")
if (HAVE_RDROBUST) {
  bws <- panel |> group_by(boundary_id) |> group_modify(function(g, key) {
    f <- rdrobust::rdrobust(y = g$n, x = g$bin, c = 0, p = 1, kernel = "triangular",
                            bwselect = "mserd", vce = "nn", cluster = g$date,
                            masspoints = "off")
    tibble::tibble(h_mserd = f$bws["h", 1], b_mserd = f$bws["b", 1])
  }) |> ungroup()
  write_csv(bws, bw_file)
} else {
  stopifnot(file.exists(bw_file))
  bws <- read_csv(bw_file, show_col_types = FALSE) |> select(boundary_id, h_mserd)
}

boundaries <- sort(unique(panel$boundary_id))
est <- lapply(boundaries, function(b) {
  g <- panel[panel$boundary_id == b, ]
  h <- bws$h_mserd[match(b, bws$boundary_id)]
  f <- rd_local_linear(g$n, g$bin, h, as.character(g$date))
  if (is.null(f)) return(NULL)
  tibble::tibble(boundary_id = b, tau = f$tau, se = f$se,
                 z = f$tau / f$se, pval = 2 * pnorm(-abs(f$tau / f$se)),
                 ci_lower = f$tau - 1.96 * f$se, ci_upper = f$tau + 1.96 * f$se,
                 bandwidth = h, n_cells = nrow(g), n_eff = f$n_eff,
                 n_eff_left = f$n_left, n_eff_right = f$n_right,
                 n_incidents = sum(g$n), n_dates = n_distinct(g$date))
}) |> bind_rows() |> arrange(tau)

## -----------------------------------------------------------------------------
## 3. Aggregation, heterogeneity, equivalence
## -----------------------------------------------------------------------------
ivw <- function(tau, se) {
  w <- 1 / se^2; t <- sum(w * tau) / sum(w); s <- 1 / sqrt(sum(w))
  list(tau = t, se = s, z = t / s, pval = 2 * pnorm(-abs(t / s)),
       lo = t - 1.96 * s, hi = t + 1.96 * s, k = length(tau))
}
P  <- ivw(est$tau, est$se)
w  <- 1 / est$se^2
Q  <- sum(w * (est$tau - P$tau)^2); dfQ <- nrow(est) - 1
I2 <- max(0, (Q - dfQ) / Q) * 100
simple_tau <- mean(est$tau); simple_se <- sqrt(sum(est$se^2)) / nrow(est)

equiv <- function(tau, se, baseline, margin = EQUIV_MARGIN) {
  d <- margin * baseline
  list(mde = 2.80 * se, mde_pct = 100 * 2.80 * se / baseline,
       upper_pct = 100 * (tau + 1.96 * se) / baseline,
       lower_pct = 100 * (tau - 1.96 * se) / baseline,
       tost_p = max(pnorm((tau - d) / se, lower.tail = TRUE),
                    pnorm((tau + d) / se, lower.tail = FALSE)),
       tightest_pct = 100 * (abs(tau) + qnorm(0.95) * se) / baseline)
}
E <- equiv(P$tau, P$se, BASELINE)

## Date-block bootstrap: resample whole dates, preserving within-date dependence.
set.seed(SEED)
dates <- sort(unique(panel$date))
idx <- split(seq_len(nrow(panel)), paste(panel$boundary_id, panel$date))
boot <- replicate(N_BOOT, {
  s  <- sample(dates, length(dates), replace = TRUE)
  tt <- vapply(boundaries, function(b) {
    keys <- paste(b, s); rows <- unlist(idx[keys[keys %in% names(idx)]], use.names = FALSE)
    if (length(rows) < 40) return(NA_real_)
    g <- panel[rows, ]
    f <- rd_local_linear(g$n, g$bin, bws$h_mserd[match(b, bws$boundary_id)], as.character(g$date))
    if (is.null(f)) NA_real_ else f$tau
  }, numeric(1))
  ok <- !is.na(tt); ww <- (1 / est$se[match(boundaries[ok], est$boundary_id)]^2)
  sum(ww * tt[ok]) / sum(ww)
})
boot <- boot[is.finite(boot)]
SE_BOOT <- sd(boot); CI_BOOT <- unname(quantile(boot, c(.025, .975)))
E_BOOT  <- equiv(P$tau, SE_BOOT, BASELINE)

## Leave-one-boundary-out
loo <- vapply(seq_len(nrow(est)), function(i) ivw(est$tau[-i], est$se[-i])$tau, numeric(1))

## -----------------------------------------------------------------------------
## 4. Robustness (Table 3)
## -----------------------------------------------------------------------------
refit <- function(p, hmap = NULL) {
  r <- lapply(sort(unique(p$boundary_id)), function(b) {
    g <- p[p$boundary_id == b, ]
    h <- if (is.null(hmap)) bws$h_mserd[match(b, bws$boundary_id)] else hmap[[b]]
    f <- rd_local_linear(g$n, g$bin, h, as.character(g$date))
    if (is.null(f)) NULL else tibble::tibble(boundary_id = b, tau = f$tau, se = f$se,
                                             n_incidents = sum(g$n))
  }) |> bind_rows()
  c(ivw(r$tau, r$se), list(n = sum(r$n_incidents), kk = nrow(r)))
}
big_ids <- est$boundary_id[est$n_incidents >= MIN_BND_N]
r_big   <- refit(panel[panel$boundary_id %in% big_ids, ])
r_donut <- refit(panel[abs(panel$bin) > DONUT_R, ])
panel100 <- inc |>
  mutate(bin = sign(running_var) * (floor(abs(running_var) / 100) * 100 + 50)) |>
  count(boundary_id, date, bin, name = "n")
grid100 <- c(seq(-MAX_D + 50, -50, by = 100), seq(50, MAX_D - 50, by = 100))
panel100 <- panel100 |> distinct(boundary_id, date) |> tidyr::crossing(bin = grid100) |>
  left_join(panel100, by = c("boundary_id","date","bin")) |> mutate(n = tidyr::replace_na(n, 0L))
r_100 <- refit(panel100)

## -----------------------------------------------------------------------------
## 5. Unit of analysis (Table 2), all rows on the same screened, zero-filled sample
## -----------------------------------------------------------------------------
pooled_ll <- function(y, x, cl, h) {
  f <- rd_local_linear(y, x, h, cl); c(f$tau, f$se, 2 * pnorm(-abs(f$tau / f$se)), f$n_eff)
}
h_pool <- median(bws$h_mserd)
u1 <- panel |> group_by(date, bin) |> summarise(n = sum(n), .groups = "drop")
u3 <- panel |> group_by(boundary_id, bin) |> summarise(n = mean(n), .groups = "drop") |>
  mutate(date = as.Date("2024-04-11"))
units <- bind_rows(
  tibble::tibble(unit = "Boundary $\\times$ date $\\times$ bin (baseline)", cells = nrow(panel),
                 tau = P$tau, se = P$se, pval = P$pval),
  {v <- pooled_ll(u1$n, u1$bin, as.character(u1$date), h_pool)
   tibble::tibble(unit = "Date $\\times$ bin, pooled over boundaries", cells = nrow(u1),
                  tau = v[1], se = v[2], pval = v[3])},
  {v <- pooled_ll(u3$n, u3$bin, as.character(u3$boundary_id), h_pool)
   tibble::tibble(unit = "Boundary $\\times$ bin, averaged over dates", cells = nrow(u3),
                  tau = v[1], se = v[2], pval = v[3])},
  ## A fourth level -- collapsing over boundaries AND dates -- leaves a single
  ## cluster, so a cluster-robust variance is not defined. It is omitted rather
  ## than reported with an infinite standard error.
  )

## -----------------------------------------------------------------------------
## 6. LaTeX tables
## -----------------------------------------------------------------------------
fmt  <- function(x, d = 3) formatC(x, format = "f", digits = d)
comma<- function(x) formatC(x, format = "d", big.mark = ",")
star <- function(p) ifelse(p < 0.01, "\\sym{***}", ifelse(p < 0.05, "\\sym{**}",
                    ifelse(p < 0.10, "\\sym{*}", "")))

t1 <- c("% Generated by 10_final_submission_estimates.R -- do not edit by hand.",
"\\begin{table}[htbp]", "\\centering",
"\\caption{Boundary-specific treatment effects of water rationing on domestic violence}",
"\\label{tab:main_results}", "\\begin{threeparttable}", "\\small",
"\\setlength{\\tabcolsep}{3.5pt}",
"\\begin{tabular}{lcccccc}", "\\toprule",
"Boundary & Effect & Std.\\ error & $p$-value & 95\\% CI & Cells & $N^{\\mathrm{eff}}$ \\\\",
"(polygon pair) & (1) & (2) & (3) & (4) & (5) & (6) \\\\", "\\midrule",
"\\addlinespace[0.1cm]",
"\\multicolumn{7}{l}{\\textit{Panel A: Individual boundaries (ordered by effect size)}} \\\\",
"\\addlinespace[0.1cm]",
sprintf("%s & $%s$%s & (%s) & %s & [$%s$, $%s$] & %s & %s \\\\",
        gsub("_", "--", est$boundary_id), fmt(est$tau, 4), star(est$pval), fmt(est$se, 4),
        fmt(est$pval, 3), fmt(est$ci_lower, 4), fmt(est$ci_upper, 4),
        comma(est$n_cells), comma(est$n_eff)),
"\\addlinespace[0.2cm]", "\\midrule", "\\addlinespace[0.1cm]",
"\\multicolumn{7}{l}{\\textit{Panel B: Pooled estimates}} \\\\", "\\addlinespace[0.1cm]",
sprintf("Inverse-variance weighted & $%s$ & (%s) & %s & [$%s$, $%s$] & %s & %s \\\\",
        fmt(P$tau,4), fmt(P$se,4), fmt(P$pval,3), fmt(P$lo,4), fmt(P$hi,4),
        comma(N_CELLS), comma(sum(est$n_eff))),
sprintf("\\quad date block bootstrap SE & & (%s) & & [$%s$, $%s$] & & \\\\",
        fmt(SE_BOOT,4), fmt(CI_BOOT[1],4), fmt(CI_BOOT[2],4)),
sprintf("Unweighted average & $%s$ & (%s) & %s & & %s & \\\\",
        fmt(simple_tau,4), fmt(simple_se,4), fmt(2*pnorm(-abs(simple_tau/simple_se)),3), comma(N_CELLS)),
"\\addlinespace[0.1cm]", "\\midrule", "\\addlinespace[0.1cm]",
"\\multicolumn{7}{l}{\\textit{Panel C: Heterogeneity}} \\\\", "\\addlinespace[0.1cm]",
sprintf("\\multicolumn{7}{l}{Cochran's $Q$ = %s ($p=%s$, %d d.f.); $I^2 = %s\\%%$} \\\\",
        fmt(Q,2), fmt(1-pchisq(Q,dfQ),3), dfQ, fmt(I2,1)),
"\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}[flushleft]", "\\small",
sprintf("\\item \\textit{Notes:} Local linear spatial regression discontinuity estimates on the balanced boundary $\\times$ date $\\times$ 50-metre-bin panel, %s cells within $\\pm$2{,}000 m of an active boundary, of which %s\\%% record no incident. The %s incidents in the deduplicated estimation sample are aggregated into these cells. A positive effect means rationing \\emph{raises} reported violence: $\\tau_b$ is the treated-side limit minus the control-side limit, and treated observations lie on the negative side of the running variable. Triangular kernel, MSE-optimal bandwidth, standard errors cluster-robust by date. Column (5) reports estimation cells, column (6) cells with positive kernel weight. The control-side mean is %s incidents per cell. \\sym{***}~$p<0.01$, \\sym{**}~$p<0.05$, \\sym{*}~$p<0.10$.",
        comma(N_CELLS), fmt(PCT_ZERO,1), comma(N_INCIDENTS), fmt(BASELINE,4)),
"\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
writeLines(t1, file.path(TAB, "table1.tex"))

t2 <- c("% Generated by 10_final_submission_estimates.R -- do not edit by hand.",
"\\begin{table}[htbp]", "\\centering",
"\\caption{Sensitivity to the unit of spatial aggregation}", "\\label{tab:units}",
"\\begin{threeparttable}", "\\small", "\\begin{tabular}{lcccc}", "\\toprule",
"Unit of analysis & Cells & Coef. & Std.\\ err. & $p$-value \\\\",
" & & (1) & (2) & (3) \\\\", "\\midrule",
sprintf("%s & %s & $%s$ & (%s) & %s \\\\", units$unit, comma(units$cells),
        fmt(units$tau,4), fmt(units$se,4), fmt(units$pval,3)),
"\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}[flushleft]", "\\small",
sprintf("\\item \\textit{Notes:} Every row uses the same balanced, zero-filled sample within $\\pm$2{,}000 m and the same local linear estimator; only the level of aggregation changes. Row 1 is the boundary-specific estimator of Table~\\ref{tab:main_results} aggregated by inverse-variance weighting; rows 2--4 pool before estimating, at a common bandwidth of %s m. Standard errors cluster-robust by date, or by boundary where dates are collapsed. Because the sample is held fixed, differences across rows reflect aggregation alone.", fmt(h_pool,0)),
"\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
writeLines(t2, file.path(TAB, "table2.tex"))

t3 <- c("% Generated by 10_final_submission_estimates.R -- do not edit by hand.",
"\\begin{table}[htbp]", "\\centering",
"\\caption{Robustness to alternative specifications}", "\\label{tab:robustness}",
"\\begin{threeparttable}", "\\small", "\\setlength{\\tabcolsep}{4.5pt}",
"\\begin{tabular}{lcccc}", "\\toprule",
" & Baseline & Large & Donut & Alternative \\\\",
" &          & boundaries & RD & bins (100 m) \\\\",
" & (1) & (2) & (3) & (4) \\\\", "\\midrule", "\\addlinespace[0.1cm]",
sprintf("Treatment effect & $%s$ & $%s$ & $%s$ & $%s$ \\\\",
        fmt(P$tau,4), fmt(r_big$tau,4), fmt(r_donut$tau,4), fmt(r_100$tau,4)),
sprintf(" & (%s) & (%s) & (%s) & (%s) \\\\",
        fmt(P$se,4), fmt(r_big$se,4), fmt(r_donut$se,4), fmt(r_100$se,4)),
"\\addlinespace[0.1cm]",
sprintf("$p$-value & %s & %s & %s & %s \\\\",
        fmt(P$pval,3), fmt(r_big$pval,3), fmt(r_donut$pval,3), fmt(r_100$pval,3)),
sprintf("95\\%% CI & [$%s$, $%s$] & [$%s$, $%s$] & [$%s$, $%s$] & [$%s$, $%s$] \\\\",
        fmt(P$lo,4), fmt(P$hi,4), fmt(r_big$lo,4), fmt(r_big$hi,4),
        fmt(r_donut$lo,4), fmt(r_donut$hi,4), fmt(r_100$lo,4), fmt(r_100$hi,4)),
"\\addlinespace[0.2cm]", "\\midrule", "\\addlinespace[0.1cm]",
sprintf("Incidents & %s & %s & %s & %s \\\\", comma(N_INCIDENTS), comma(r_big$n),
        comma(r_donut$n), comma(r_100$n)),
sprintf("Boundaries & %d & %d & %d & %d \\\\", P$k, r_big$kk, r_donut$kk, r_100$kk),
"Bin width (m) & 50 & 50 & 50 & 100 \\\\",
sprintf("Exclusion window & -- & -- & $\\pm$%d m & -- \\\\", DONUT_R),
sprintf("Sample restriction & none & $N_b \\geq %d$ & none & none \\\\", MIN_BND_N),
"\\addlinespace[0.2cm]", "\\midrule", "\\addlinespace[0.1cm]",
sprintf("\\multicolumn{5}{l}{Leave-one-boundary-out range of $\\hat\\tau_{\\mathrm{IV}}$: $[%s,\\,%s]$} \\\\",
        fmt(min(loo),4), fmt(max(loo),4)),
"\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}[flushleft]", "\\small",
"\\item \\textit{Notes:} Pooled inverse-variance weighted treatment effects on the balanced, zero-filled panel. Column (2) drops boundaries with fewer than 500 incidents within 2{,}000 m. Column (3) excludes cells within 100 m of the boundary. Column (4) rebuilds the panel on 100-metre bins. All columns use a triangular kernel, MSE-optimal bandwidths and date-clustered standard errors. The final row reports the range of the pooled estimate across the leave-one-boundary-out recomputations.",
"\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
writeLines(t3, file.path(TAB, "table3.tex"))

srow <- function(lab, v, d = 0) sprintf("%s & %s & %s & %s & %s & %s \\\\", lab,
  formatC(mean(v), format="f", digits=d, big.mark=","), formatC(sd(v), format="f", digits=d, big.mark=","),
  formatC(min(v),  format="f", digits=d, big.mark=","), formatC(median(v), format="f", digits=d, big.mark=","),
  formatC(max(v),  format="f", digits=d, big.mark=","))
tc <- c("% Generated by 10_final_submission_estimates.R -- do not edit by hand.",
"\\begin{table}[htbp]", "\\centering",
"\\caption{Summary statistics: boundary-level characteristics}", "\\label{tab:summary_stats}",
"\\begin{threeparttable}", "\\small", "\\begin{tabular}{lcccccc}", "\\toprule",
"Variable & Mean & Std.\\ dev. & Min & Median & Max \\\\", "\\midrule", "\\addlinespace[0.1cm]",
"\\multicolumn{6}{l}{\\textit{Panel A: Sample composition}} \\\\", "\\addlinespace[0.1cm]",
srow("Incidents within 2,000 m of the boundary", est$n_incidents),
srow("Dates on which the boundary is active", est$n_dates),
srow("Estimation cells (date $\\times$ bin)", est$n_cells),
srow("Cells with positive kernel weight", est$n_eff),
"\\addlinespace[0.2cm]", "\\midrule", "\\addlinespace[0.1cm]",
"\\multicolumn{6}{l}{\\textit{Panel B: Bandwidth}} \\\\", "\\addlinespace[0.1cm]",
srow("MSE-optimal bandwidth (m)", est$bandwidth),
"\\addlinespace[0.2cm]", "\\midrule", "\\addlinespace[0.1cm]",
"\\multicolumn{6}{l}{\\textit{Panel C: Boundary-specific treatment effects}} \\\\", "\\addlinespace[0.1cm]",
srow("Point estimate ($\\hat{\\tau}_b$)", est$tau, 4),
srow("Standard error", est$se, 4),
srow("Absolute effect size", abs(est$tau), 4),
"\\bottomrule", "\\end{tabular}", "\\begin{tablenotes}[flushleft]", "\\small",
sprintf("\\item \\textit{Notes:} Statistics across the %d boundaries that generate identifying variation, unweighted. Incidents are counted within 2{,}000 m of the assigned boundary and total %s. Estimation cells are boundary $\\times$ date $\\times$ 50-metre-bin units of the balanced zero-filled panel and total %s; %s of them carry positive kernel weight at the MSE-optimal bandwidth. A boundary is active only on dates when one of its two polygons is rationed, which is why the date counts are far below the %d-day study window.",
        nrow(est), comma(N_INCIDENTS), comma(N_CELLS), comma(sum(est$n_eff)), 325L),
"\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
writeLines(tc, file.path(TAB, "table_c1.tex"))

## -----------------------------------------------------------------------------
## 7. Machine-readable outputs and a manifest of every in-text statistic
## -----------------------------------------------------------------------------
write_csv(est, file.path(OUT, "boundary_estimates_final.csv"))
stats <- tibble::tribble(
  ~quantity, ~value,
  "n_incidents",                 N_INCIDENTS,
  "n_cells",                     N_CELLS,
  "pct_zero_cells",              PCT_ZERO,
  "n_boundaries",                nrow(est),
  "n_boundary_dates",            n_distinct(paste(panel$boundary_id, panel$date)),
  "control_baseline_per_cell",   BASELINE,
  "pooled_tau",                  P$tau,
  "pooled_se_meta",              P$se,
  "pooled_se_bootstrap",         SE_BOOT,
  "pooled_se_ratio",             SE_BOOT / P$se,
  "pooled_pval",                 P$pval,
  "pooled_ci_lower",             P$lo,
  "pooled_ci_upper",             P$hi,
  "boot_ci_lower",               CI_BOOT[1],
  "boot_ci_upper",               CI_BOOT[2],
  "simple_tau",                  simple_tau,
  "simple_se",                   simple_se,
  "cochran_Q",                   Q,
  "cochran_Q_pval",              1 - pchisq(Q, dfQ),
  "I2_pct",                      I2,
  "mde_meta",                    E$mde,
  "mde_pct_meta",                E$mde_pct,
  "mde_pct_boot",                E_BOOT$mde_pct,
  "upper_ci_pct_meta",           E$upper_pct,
  "upper_ci_pct_boot",           E_BOOT$upper_pct,
  "tost_p_meta",                 E$tost_p,
  "tost_p_boot",                 E_BOOT$tost_p,
  "tightest_equiv_pct_meta",     E$tightest_pct,
  "tightest_equiv_pct_boot",     E_BOOT$tightest_pct,
  "loo_min",                     min(loo),
  "loo_max",                     max(loo),
  "bandwidth_mean",              mean(est$bandwidth),
  "bandwidth_min",               min(est$bandwidth),
  "bandwidth_max",               max(est$bandwidth),
  "n_eff_total",                 sum(est$n_eff),
  "boot_B",                      length(boot))
write_csv(stats, file.path(OUT, "manuscript_statistics.csv"))

cat("\n================ FINAL ESTIMATES ================\n")
cat(sprintf("Incidents %s | cells %s (%.1f%% zero) | boundaries %d | baseline %.4f/cell\n",
            comma(N_INCIDENTS), comma(N_CELLS), PCT_ZERO, nrow(est), BASELINE))
cat(sprintf("Pooled tau (treated - control) = %+.5f\n", P$tau))
cat(sprintf("  meta-analytic SE %.5f  p = %.3f  CI [%+.5f, %+.5f]\n", P$se, P$pval, P$lo, P$hi))
cat(sprintf("  bootstrap SE     %.5f  (ratio %.3f)  percentile CI [%+.5f, %+.5f]  B = %d\n",
            SE_BOOT, SE_BOOT/P$se, CI_BOOT[1], CI_BOOT[2], length(boot)))
cat(sprintf("Q = %.2f (p = %.3f), I2 = %.1f%%\n", Q, 1-pchisq(Q,dfQ), I2))
cat(sprintf("Equivalence, meta-analytic : MDE %.1f%%  |CI| upper %+.1f%%  tightest %.1f%%  TOST p = %.4f\n",
            E$mde_pct, E$upper_pct, E$tightest_pct, E$tost_p))
cat(sprintf("Equivalence, bootstrap     : MDE %.1f%%  |CI| upper %+.1f%%  tightest %.1f%%  TOST p = %.4f\n",
            E_BOOT$mde_pct, E_BOOT$upper_pct, E_BOOT$tightest_pct, E_BOOT$tost_p))
cat(sprintf("Leave-one-out range: [%+.5f, %+.5f]\n", min(loo), max(loo)))
cat("Tables written to", TAB, "\n")
