#!/usr/bin/env Rscript
## =============================================================================
## 12_randomization_inference.R
##
## EXACT DESIGN-BASED INFERENCE for the pooled spatial RD estimate, computed
## strictly on the balanced, zero-filled boundary x date x bin panel.
##
## =============================================================================
## WHY THE CALENDAR RE-PHASING IN THE REVISION BRIEF CANNOT BE USED
## =============================================================================
## Two candidate re-phasings were implemented and both are degenerate on this
## design. Both diagnostics are retained here because they are reported in the
## manuscript appendix and in REVISION_NOTES.md.
##
## (i) RE-PHASING THE CITYWIDE CALENDAR (the literal request: 324 non-identity
##     rotations of the 325-day window). The rotation is a deterministic 9-day
##     cycle over 9 shifts. A boundary between polygons in shifts s_p and s_q is
##     ACTIVE on a date only when exactly one of the two is rationed. Shifting the
##     calendar by k days maps every treated shift s to s + k (mod 9), so such a
##     boundary is active only when k = 0 or k = (s_q - s_p) mod 9 -- two phases in
##     nine. The diagnostic on this sample: the realised calendar gives a typical
##     boundary 70 active dates, the mean across the 324 shifted calendars is 14.1,
##     and 49.1% of boundary-phase cells have NO active dates at all. The observed
##     statistic would be compared against a distribution computed from a small,
##     data-starved fraction of the design, whose dispersion (SD 0.0362) reflects
##     noise in tiny subsamples rather than the sampling variability of the observed
##     estimator (bootstrap SE 0.0081). The RI figure previously reported by
##     06_design_checks.py is subject to this defect and is withdrawn.
##
## (ii) RE-PHASING THE SIDE SEQUENCE WITHIN BOUNDARY (composition preserving).
##     Because the two polygons of a boundary sit in two fixed shifts of a regular
##     9-day rotation, the realised treated-side sequence is an exact alternation
##     p, q, p, q, ... broken only where the calendar itself was interrupted. A
##     cyclic shift of an alternating sequence returns either the sequence itself
##     (even shifts) or its complement (odd shifts), so the reference distribution
##     collapses onto {+tau, -tau} and has no power.
##
## =============================================================================
## THE RANDOMIZATION ACTUALLY USED: SIGN-FLIP RELABELING
## =============================================================================
## The randomness the design can credibly invoke is not WHEN a shift is rationed --
## that is a published deterministic calendar -- but WHICH member of an adjacent
## polygon pair was placed in the earlier of the two shifts when the city drew up
## its hydraulic sectorisation. Relabeling the two shifts of a boundary exchanges
## the treated and control sides on every active date of that boundary, and by the
## symmetry of the bin grid and the common bandwidth it maps tau_b to -tau_b
## exactly. The panel, the incidents, the distances, the bins, the active dates and
## the treated/control day counts are all preserved.
##
##  DESIGN A (primary, exact by complete enumeration)
##    eps_b in {+1,-1} independently for each of the 15 boundaries: all
##    2^15 = 32,768 relabelings, 32,767 of them non-identity. No simulation, no
##    asymptotics and no model of the spatial covariance is involved.
##
##  DESIGN B (secondary, finer, exact under uniform sampling of the group)
##    The 9-day rotation partitions the window into complete cycles; within a cycle
##    each boundary is rationed once on each side. Flipping independently by
##    boundary AND cycle randomises which side of the pair is rationed first within
##    each cycle -- the local-randomization RD assumption stated at the cycle level.
##    The group has 2^(sum_b n_cycles_b) elements and is sampled uniformly.
##
## Both statistics are computed in closed form. Writing a_b for the row of
## (n_b G_b)^{-1} U_b that returns the discontinuity coefficient, the estimator is
## linear in the cell counts, so for any flip pattern f over a boundary's blocks
##      tau_b(f) = tau_b(0) + f' v_b ,   v_b = (A0_b - A1_b) (aN_b - aP_b),
## where A1_b and A0_b are the block-aggregated count matrices on the realised
## treated and control sides. The whole reference distribution is therefore a
## matrix product.
##
## =============================================================================
## CONFIDENCE INTERVAL BY INVERSION
## =============================================================================
## Linearity also gives the interval in closed form: under the sharp null of a
## constant effect tau0 the adjusted outcome is y - tau0*D_obs and the relabeled
## statistic is exactly T_g(y) - tau0*T_g(D_obs). Computing T_g(y) and T_g(D_obs)
## once per relabeling delivers the entire inverted interval. The 95% randomization
## interval is the set of tau0 not rejected at the 5% level; its outer edges are
## the exact RI equivalence bounds, which require no normal approximation and no
## 2.80 x SE minimum-detectable-effect rule.
##
## =============================================================================
## WHAT IS AND IS NOT EXACT
## =============================================================================
## The test is exact CONDITIONAL on the assignment of incidents to boundaries. That
## assignment was made under the realised calendar (each incident is attached to
## its nearest ACTIVE boundary), so a fully unconditional test would re-run the
## assignment under every relabeling, which requires the pre-screen geocoded
## extract. This is stated in the manuscript.
## =============================================================================

suppressPackageStartupMessages({library(dplyr); library(readr)})
ROOT <- Sys.getenv("ROOT", "."); DER <- file.path(ROOT,"data","derived")
PUB  <- file.path(ROOT,"data","public"); OUT <- file.path(ROOT,"output")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)

BIN_W <- 50; MAX_D <- 2000
MAGS  <- seq(BIN_W/2, MAX_D-BIN_W/2, by=BIN_W)     # 40 bin midpoints
CYCLE <- 9L                                        # days in the rationing rotation
N_MC  <- 200000L                                   # uniform draws for design B
SEED  <- 20240411L

## ---- inputs -----------------------------------------------------------------
inc <- read_csv(file.path(DER,"rd_estimation_sample.csv"), show_col_types=FALSE) |>
  filter(in_estimation_sample == 1, is_duplicate_record == 0) |>
  mutate(date = as.Date(date),
         mag  = floor(abs(running_var)/BIN_W)*BIN_W + BIN_W/2)
panel <- read_csv(file.path(DER,"panel_balanced_zerofilled.csv"), show_col_types=FALSE) |>
  mutate(date = as.Date(date))
est <- read_csv(file.path(OUT,"boundary_estimates_final.csv"), show_col_types=FALSE)

stopifnot(nrow(inc) == sum(panel$n))               # 20,527 incidents, both sources

bnds <- sort(unique(inc$boundary_id)); NB <- length(bnds)
Hb   <- setNames(est$bandwidth[match(bnds, est$boundary_id)], bnds)
Wb   <- setNames(1/est$se[match(bnds, est$boundary_id)]^2,   bnds)
Tb   <- setNames(est$tau[match(bnds, est$boundary_id)],      bnds)
DAY0 <- min(panel$date)

## ---- the linear functional a_b that extracts the discontinuity ---------------
grid_x <- c(-rev(MAGS), MAGS)                      # 80 signed bin midpoints
avec <- function(h, n) {
  w <- pmax(0, 1 - abs(grid_x)/h); D <- as.numeric(grid_x < 0)
  Z <- cbind(1, D, grid_x, grid_x*D)               # 80 x 4
  U <- t(Z * w); G <- crossprod(Z * sqrt(w))       # U: 4 x 80 ; G: 4 x 4 per date
  as.numeric(solve(n*G, U)[2, ])                   # length 80
}

## ---- per boundary: realised treated / control side count matrices ------------
pre <- lapply(bnds, function(b) {
  pq <- as.integer(strsplit(b, "_")[[1]]); p <- pq[1]
  d  <- inc[inc$boundary_id == b, ]
  dts <- sort(unique(panel$date[panel$boundary_id == b]))   # active dates from panel
  n  <- length(dts)
  C1 <- matrix(0, n, length(MAGS)); C0 <- C1            # treated side / control side
  di <- match(d$date, dts); mi <- match(d$mag, MAGS)
  stopifnot(!anyNA(di), !anyNA(mi))
  tr <- d$treated == 1
  for (j in which(tr))  C1[di[j], mi[j]] <- C1[di[j], mi[j]] + 1
  for (j in which(!tr)) C0[di[j], mi[j]] <- C0[di[j], mi[j]] + 1
  a  <- avec(Hb[[b]], n); aN <- rev(a[1:40]); aP <- a[41:80]
  cyc <- as.integer((as.integer(dts - DAY0)) %/% CYCLE)
  ucy <- sort(unique(cyc)); ci <- match(cyc, ucy)
  A1 <- rowsum(C1, ci); A0 <- rowsum(C0, ci)            # block-aggregated counts
  m  <- as.numeric(table(ci))                           # active dates per cycle
  list(n=n, aN=aN, aP=aP, ncyc=length(ucy), m=m,
       tau0 = sum(aN*colSums(C1)) + sum(aP*colSums(C0)),
       v    = as.numeric((A0 - A1) %*% (aN - aP)),      # length ncyc
       sN   = sum(aN), sP = sum(aP))
})
names(pre) <- bnds

## ---- verification against the balanced panel and against Table 1 -------------
recon <- do.call(rbind, lapply(bnds, function(b) {
  d <- inc[inc$boundary_id == b, ]
  dts <- sort(unique(panel$date[panel$boundary_id == b]))
  C1 <- matrix(0, length(dts), length(MAGS)); C0 <- C1
  di <- match(d$date, dts); mi <- match(d$mag, MAGS); tr <- d$treated == 1
  for (j in which(tr))  C1[di[j], mi[j]] <- C1[di[j], mi[j]] + 1
  for (j in which(!tr)) C0[di[j], mi[j]] <- C0[di[j], mi[j]] + 1
  data.frame(boundary_id=b,
             date=rep(dts, times=2*length(MAGS)),
             bin =rep(c(-rev(MAGS), MAGS), each=length(dts)),
             nc  =c(as.vector(C1[, rev(seq_along(MAGS))]), as.vector(C0)))
}))
cmp <- panel |> inner_join(recon, by=c("boundary_id","date","bin"))
stopifnot(nrow(cmp) == nrow(panel), all(cmp$n == cmp$nc))
cat(sprintf("panel reconstruction verified : %s cells, %s incidents, identical\n",
            format(nrow(panel), big.mark=","), format(sum(panel$n), big.mark=",")))

tau_obs_b <- vapply(pre, function(P) P$tau0, numeric(1))
cat(sprintf("boundary tau vs Table 1 file  : max abs difference %.3e\n",
            max(abs(tau_obs_b - Tb))))
stopifnot(max(abs(tau_obs_b - Tb)) < 1e-9)
## a full flip must send tau_b to exactly -tau_b
flip_all <- vapply(bnds, function(b) pre[[b]]$tau0 + sum(pre[[b]]$v), numeric(1))
cat(sprintf("full-relabel antisymmetry     : max |tau(flip) + tau| = %.3e\n",
            max(abs(flip_all + tau_obs_b))))
stopifnot(max(abs(flip_all + tau_obs_b)) < 1e-9)

WSUM  <- sum(Wb)
obs_y <- sum(Wb*tau_obs_b)/WSUM
## the statistic applied to the realised treatment indicator
tauD_b <- vapply(bnds, function(b) pre[[b]]$n * pre[[b]]$sN, numeric(1))
obs_d  <- sum(Wb*tauD_b)/WSUM
cat(sprintf("observed pooled tau           : %+.8f\n", obs_y))
cat(sprintf("T(D) at the identity          : %.10f\n", obs_d))
stopifnot(abs(obs_d - 1) < 1e-8)

## =============================================================================
## DESIGN A -- complete enumeration of the 2^15 boundary relabelings
## =============================================================================
## eps_b = -1 flips boundary b entirely: tau_b -> -tau_b and tauD_b -> -tauD_b.
E <- matrix(0L, 2^NB, NB)
for (j in seq_len(NB)) E[, j] <- rep(c(0L,1L), each = 2^(j-1), length.out = 2^NB)
S <- 1 - 2*E                                        # +1 keep, -1 flip
Ay <- as.numeric(S %*% (Wb*tau_obs_b)) / WSUM
Ad <- as.numeric(S %*% (Wb*tauD_b))   / WSUM
stopifnot(abs(Ay[1] - obs_y) < 1e-12, abs(Ad[1] - obs_d) < 1e-12)
Ay <- Ay[-1]; Ad <- Ad[-1]                          # drop the identity

## =============================================================================
## DESIGN B -- uniform draws from the boundary x cycle relabeling group
## =============================================================================
ncy <- vapply(pre, function(P) P$ncyc, integer(1))
set.seed(SEED)
acc_y <- matrix(0, N_MC - 1L, NB); acc_d <- matrix(0, N_MC - 1L, NB)
for (i in seq_len(NB)) {
  P <- pre[[i]]
  Fm <- matrix(rbinom((N_MC-1L)*P$ncyc, 1L, 0.5), N_MC-1L, P$ncyc)
  acc_y[, i] <- P$tau0 + as.numeric(Fm %*% P$v)
  K <- as.numeric(Fm %*% P$m)                       # flipped active dates
  acc_d[, i] <- P$n*P$sN + K*(P$sP - P$sN)
}
By <- as.numeric(acc_y %*% Wb)/WSUM
Bd <- as.numeric(acc_d %*% Wb)/WSUM

## =============================================================================
## EXACT p-VALUES AND INVERTED INTERVALS
## =============================================================================
p_two <- function(null, o) (1 + sum(abs(null) >= abs(o))) / (1 + length(null))
pA <- p_two(Ay, obs_y); pB <- p_two(By, obs_y)
pctA <- 100*mean(Ay <= obs_y); pctB <- 100*mean(By <= obs_y)

invert <- function(ny, nd, alpha = 0.05) {
  span <- max(12*sd(ny), 40*abs(obs_y))
  grid <- seq(obs_y - span, obs_y + span, length.out = 20001)
  keep <- vapply(grid, function(t0) {
    s_o <- obs_y - t0*obs_d; s_n <- ny - t0*nd
    (1 + sum(abs(s_n) >= abs(s_o))) / (1 + length(s_n)) >= alpha
  }, logical(1))
  stopifnot(any(keep), !keep[1], !keep[length(keep)])
  range(grid[keep])
}
one_sided_up <- function(ny, nd, alpha = 0.05) {
  span <- max(12*sd(ny), 40*abs(obs_y))
  grid <- seq(obs_y - span, obs_y + span, length.out = 20001)
  keep <- vapply(grid, function(t0) {
    s_o <- obs_y - t0*obs_d; s_n <- ny - t0*nd
    (1 + sum(s_n <= s_o)) / (1 + length(s_n)) >= alpha
  }, logical(1))
  max(grid[keep])
}
ciA <- invert(Ay, Ad); ciB <- invert(By, Bd)
upA <- one_sided_up(Ay, Ad); upB <- one_sided_up(By, Bd)

## =============================================================================
## BIN-WIDTH-INVARIANT POLICY MAGNITUDES
## =============================================================================
## tau is a discontinuity in expected incidents per 50 m bin per boundary-day. The
## bin width cancels once the estimate is aggregated over the band the estimator
## actually uses: multiplying by the number of bins inside the bandwidth converts
## it into incidents per boundary-day, which is invariant to BIN_W because halving
## the bins halves tau and doubles their number. Multiplying again by the number of
## boundaries active on a rationing day gives the implied citywide daily total.
BASE  <- mean(panel$n[panel$bin > 0])                  # control-side mean per cell
bd    <- panel |> distinct(boundary_id, date)
n_bd  <- nrow(bd); n_dates <- n_distinct(bd$date)
bnd_per_day <- n_bd / n_dates
inc_per_day_perimeter <- sum(panel$n) / n_dates
CITYWIDE_DV_PER_DAY <- as.numeric(Sys.getenv("CITYWIDE_DV_PER_DAY", "302.8"))

bins_in_bw <- vapply(bnds, function(b) sum(MAGS < Hb[[b]]), numeric(1))
BINS_W <- sum(Wb*bins_in_bw)/WSUM
band_m <- BINS_W * BIN_W

mag <- function(t) c(per_bin_day   = t,
                     per_bnd_day   = t*BINS_W,
                     citywide_day  = t*BINS_W*bnd_per_day,
                     pct_perimeter = 100*t*BINS_W*bnd_per_day/inc_per_day_perimeter,
                     pct_citywide  = 100*t*BINS_W*bnd_per_day/CITYWIDE_DV_PER_DAY,
                     pct_baseline  = 100*t/BASE)
## PRIMARY magnitudes use design B (the finer, cycle-level relabeling), which is
## the randomization the local-randomization RD framing actually implies; design A
## is reported alongside as the conservative complete-enumeration bracket.
m_obs <- mag(obs_y)
m_up  <- mag(ciB[2]); m_lo  <- mag(ciB[1])      # primary
c_up  <- mag(ciA[2]); c_lo  <- mag(ciA[1])      # conservative

## ---- heterogeneity on the balanced panel ------------------------------------
Qstat <- sum(Wb*(tau_obs_b - obs_y)^2); dfQ <- NB - 1
pQ <- pchisq(Qstat, dfQ, lower.tail=FALSE); I2 <- max(0, 100*(Qstat - dfQ)/Qstat)
## Design-based counterpart. Under design A, sum_b w_b tau_b^2 is invariant to the
## relabeling, so Q(eps) = sum_b w_b tau_b^2 - W * taubar(eps)^2 is a strictly
## DECREASING function of |taubar|: the boundary-relabeling null cannot separate
## "excess dispersion" from "a pooled mean near zero" -- they are the same
## statistic, and the observed Q is large precisely BECAUSE the pooled estimate is
## almost exactly zero. Q is therefore tested against the finer design B, under
## which the individual tau_b genuinely move.
Qb   <- rowSums(sweep(acc_y - as.vector(By), 2, Wb, `*`) * (acc_y - as.vector(By)))
pQ_ri <- (1 + sum(Qb >= Qstat)) / (1 + length(Qb))
I2_ri <- 100*max(0, 1 - median(Qb)/Qstat)

## Because the design-based Q test rejects at the 5% level, the boundary-specific
## estimates are tested individually against their own exact null, and a max-|z|
## statistic gives a familywise-exact test that controls for the 15 comparisons
## without a Bonferroni penalty (Westfall-Young, on the same draws).
se_b   <- est$se[match(bnds, est$boundary_id)]
Zobs   <- tau_obs_b/se_b
Znull  <- sweep(acc_y, 2, se_b, `/`)
p_b    <- (1 + colSums(abs(Znull) >= matrix(abs(Zobs), nrow(Znull), NB, byrow=TRUE))) /
          (1 + nrow(Znull))
p_maxT <- (1 + sum(apply(abs(Znull), 1, max) >= max(abs(Zobs)))) / (1 + nrow(Znull))
bnd_tab <- tibble::tibble(boundary_id = bnds, tau = tau_obs_b, se = se_b,
                          z = Zobs, p_ri = as.numeric(p_b),
                          n_incidents = est$n_incidents[match(bnds, est$boundary_id)])
write_csv(bnd_tab, file.path(OUT,"ri_boundary_pvalues.csv"))
rm(acc_y, acc_d, Znull); invisible(gc())

## ---- write ------------------------------------------------------------------
out <- tibble::tibble(
  quantity = c("ri_observed_tau","ri_p_exact","ri_p_cycle","ri_percentile",
               "ri_n_relabelings","ri_n_cycle_draws",
               "ri_ci_lower","ri_ci_upper","ri_ci_lower_cycle","ri_ci_upper_cycle",
               "ri_one_sided_upper","ri_equiv_margin_abs",
               "ri_equiv_margin_pct_baseline","ri_null_sd","ri_null_sd_cycle",
               "meta_analytic_se","bootstrap_se",
               "control_baseline_per_cell","bins_in_bandwidth","band_width_m",
               "boundaries_per_rationing_day","incidents_per_day_perimeter",
               "citywide_dv_per_day",
               "effect_per_boundary_day","effect_citywide_per_day",
               "effect_pct_of_perimeter_flow","effect_pct_of_citywide_flow",
               "ci_upper_citywide_per_day","ci_upper_pct_of_perimeter_flow",
               "ci_upper_pct_of_citywide_flow",
               "ci_lower_citywide_per_day","ci_lower_pct_of_citywide_flow",
               "equiv_citywide_per_day","equiv_pct_of_citywide_flow",
               "equiv_margin_abs_conservative","equiv_margin_pct_baseline_conservative",
               "ci_upper_citywide_per_day_conservative",
               "ci_lower_citywide_per_day_conservative",
               "equiv_citywide_per_day_conservative",
               "Q_stat","Q_df","Q_pvalue","I2_pct","Q_pvalue_ri","I2_ri_pct","p_familywise_maxT"),
  value = c(obs_y, pA, pB, pctA,
            length(Ay)+1, length(By),
            ciA[1], ciA[2], ciB[1], ciB[2],
            upA, max(abs(ciB)),
            100*max(abs(ciB))/BASE, sd(Ay), sd(By),
            1/sqrt(WSUM), 0.008130,
            BASE, BINS_W, band_m,
            bnd_per_day, inc_per_day_perimeter, CITYWIDE_DV_PER_DAY,
            m_obs[["per_bnd_day"]], m_obs[["citywide_day"]],
            m_obs[["pct_perimeter"]], m_obs[["pct_citywide"]],
            m_up[["citywide_day"]], m_up[["pct_perimeter"]], m_up[["pct_citywide"]],
            m_lo[["citywide_day"]], m_lo[["pct_citywide"]],
            mag(max(abs(ciB)))[["citywide_day"]], mag(max(abs(ciB)))[["pct_citywide"]],
            max(abs(ciA)), 100*max(abs(ciA))/BASE,
            c_up[["citywide_day"]], c_lo[["citywide_day"]],
            mag(max(abs(ciA)))[["citywide_day"]],
            Qstat, dfQ, pQ, I2, pQ_ri, I2_ri, p_maxT))
write_csv(out, file.path(OUT,"randomization_inference.csv"))
write_csv(tibble::tibble(design="A_boundary_relabeling",
                         draw=seq_along(Ay), tau=Ay, tau_D=Ad),
          file.path(OUT,"ri_null_distribution.csv"))

## =============================================================================
## FIGURE 8 -- the exact randomization null distribution
## =============================================================================
suppressPackageStartupMessages({library(ggplot2); library(tidyr)})
FIG <- file.path(ROOT,"figures"); dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
PURPLE <- "#660099"; LILAC <- "#B280CC"
fdat <- bind_rows(
  tibble::tibble(design = sprintf("(a) Boundary relabelings: all 2^15 = 32,768 (exact enumeration)"),
                 tau = Ay),
  tibble::tibble(design = sprintf("(b) Boundary x rotation-cycle relabelings: 199,999 draws"),
                 tau = By))
p8 <- ggplot(fdat, aes(tau)) +
  geom_histogram(bins = 90, fill = LILAC, colour = NA) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = obs_y, colour = PURPLE, linewidth = 0.9) +
  facet_wrap(~design, ncol = 1, scales = "free") +
  labs(x = expression(hat(tau)[IV]~"under relabeling (incidents per boundary "%*%" date "%*%" bin cell)"),
       y = "Relabelings") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.text = element_text(hjust = 0))
ggsave(file.path(FIG,"Figure_8.pdf"), p8, width = 6.5, height = 5.4, device = cairo_pdf)

## ---- report ------------------------------------------------------------------
cat("\n=============== EXACT RANDOMIZATION INFERENCE ===============\n")
cat(sprintf("balanced zero-filled panel      : %d boundary-dates x %d bins = %s cells\n",
            n_bd, 2*length(MAGS), format(nrow(panel), big.mark=",")))
cat(sprintf("observed pooled tau             : %+.6f incidents per 50 m bin per boundary-day\n", obs_y))
cat(sprintf("\n(A) boundary relabelings -- CONSERVATIVE, complete enumeration of 2^%d = %s\n",
            NB, format(2^NB, big.mark=",")))
cat(sprintf("    exact two-sided p-value     : %.4f   (%s non-identity relabelings)\n",
            pA, format(length(Ay), big.mark=",")))
cat(sprintf("    observed sits at the        : %.1fth percentile of the null\n", pctA))
cat(sprintf("    null SD                     : %.6f\n", sd(Ay)))
cat(sprintf("    95%% randomization interval  : [%+.5f, %+.5f]\n", ciA[1], ciA[2]))
cat(sprintf("    one-sided 95%% upper bound   : %+.5f\n", upA))
cat(sprintf("\n(B) boundary x rotation-cycle relabelings -- PRIMARY, %s uniform draws\n",
            format(N_MC-1L, big.mark=",")))
cat(sprintf("    exact two-sided p-value     : %.4f   (%.1fth percentile, null SD %.6f)\n",
            pB, pctB, sd(By)))
cat(sprintf("    95%% randomization interval  : [%+.5f, %+.5f]\n", ciB[1], ciB[2]))
cat(sprintf("\nfor comparison, meta-analytic SE %.6f ; date-block bootstrap SE %.6f\n",
            1/sqrt(WSUM), 0.008130))
cat(sprintf("exact RI equivalence margin     : %.5f = %.1f%% of the control-side mean (%.4f)\n",
            max(abs(ciB)), 100*max(abs(ciB))/BASE, BASE))
cat(sprintf("  conservative (design A)       : %.5f = %.1f%% of the control-side mean\n",
            max(abs(ciA)), 100*max(abs(ciA))/BASE))
cat("\n---------------- HETEROGENEITY (balanced panel) ----------------\n")
cat(sprintf("Cochran Q = %.2f on %d df, p = %.3f ; I2 = %.1f%%\n", Qstat, dfQ, pQ, I2))
cat(sprintf("design-based p for Q (cycle relabeling)         : %.4f\n", pQ_ri))
cat(sprintf("design-based I2 (1 - median null Q / observed Q) : %.1f%%\n", I2_ri))
cat(sprintf("familywise-exact max|z| p-value across 15 boundaries : %.4f\n", p_maxT))
print(as.data.frame(bnd_tab |> arrange(p_ri)), digits=3, row.names=FALSE)
cat("\n---------------- BIN-WIDTH-INVARIANT POLICY MAGNITUDES ----------------\n")
cat(sprintf("weighted bins inside bandwidth  : %.2f  (band = %.0f m each side)\n", BINS_W, band_m))
cat(sprintf("boundaries active per day       : %.2f\n", bnd_per_day))
cat(sprintf("incidents/day within perimeters : %.1f\n", inc_per_day_perimeter))
cat(sprintf("citywide DV reports per day     : %.1f\n", CITYWIDE_DV_PER_DAY))
cat(sprintf("point estimate  : %+.4f per boundary-day  ->  %+.4f citywide incidents per rationing day\n",
            m_obs[["per_bnd_day"]], m_obs[["citywide_day"]]))
cat(sprintf("                : %+.3f%% of the perimeter daily flow ; %+.4f%% of citywide reports\n",
            m_obs[["pct_perimeter"]], m_obs[["pct_citywide"]]))
cat(sprintf("95%% RI bounds   : [%+.4f, %+.4f] citywide incidents per rationing day\n",
            m_lo[["citywide_day"]], m_up[["citywide_day"]]))
cat(sprintf("                : [%+.3f%%, %+.3f%%] of citywide daily reports\n",
            m_lo[["pct_citywide"]], m_up[["pct_citywide"]]))
cat(sprintf("conservative RI bounds (A) : [%+.4f, %+.4f] citywide incidents per rationing day\n",
            c_lo[["citywide_day"]], c_up[["citywide_day"]]))
cat(sprintf("equivalence     : the data rule out effects larger than %.3f citywide incidents\n",
            mag(max(abs(ciB)))[["citywide_day"]]))
cat(sprintf("                  per rationing day, i.e. %.3f%% of citywide daily reports\n",
            mag(max(abs(ciB)))[["pct_citywide"]]))
cat(sprintf("                  (conservative: %.3f incidents, %.3f%%)\n",
            mag(max(abs(ciA)))[["citywide_day"]], mag(max(abs(ciA)))[["pct_citywide"]]))
cat("=============================================================\n")

## =============================================================================
## LATEX OUTPUT
##   (i) Table 1 Panels B and C are augmented with the exact design-based rows,
##       so that the primary inference reported in the text is the one in the
##       table. Anchors are the deterministic strings written by step 10.
##   (ii) Table 4 (design-based inference) is regenerated in full.
## =============================================================================
TAB <- file.path(ROOT, "tables")
f1  <- file.path(TAB, "table1.tex")
t1  <- readLines(f1, warn = FALSE)
## idempotent: strip any rows a previous run of this step inserted
t1  <- t1[!grepl("exact randomization test|conservative relabeling|Design-based \\$Q\\$ test", t1)]
t1  <- sub("Panel B reports, alongside.*?The control-side mean is", "The control-side mean is", t1)
fmt <- function(x, d = 4) formatC(x, format = "f", digits = d)

a_b <- grep("^\\\\quad date block bootstrap SE", t1)
a_c <- grep("Cochran's \\$Q\\$", t1)
stopifnot(length(a_b) == 1, length(a_c) == 1)

t1[a_b] <- paste0(t1[a_b], "\n",
  sprintf("\\quad exact randomization test & & & %s & [$%s$, $%s$] & & \\\\",
          fmt(pB, 3), fmt(ciB[1]), fmt(ciB[2])), "\n",
  sprintf("\\quad \\emph{conservative relabeling} & & & %s & [$%s$, $%s$] & & \\\\",
          fmt(pA, 3), fmt(ciA[1]), fmt(ciA[2])))
t1[a_c] <- paste0(t1[a_c], "\n",
  sprintf("\\multicolumn{7}{l}{Design-based $Q$ test ($p=%s$); familywise exact $\\max|z|$ ($p=%s$)} \\\\",
          fmt(pQ_ri, 3), fmt(p_maxT, 3)))

nt <- grep("The control-side mean is", t1)
stopifnot(length(nt) == 1)
t1[nt] <- sub("The control-side mean is",
  sprintf(paste0("Panel B reports, alongside the model-based standard errors, the exact ",
                 "design-based test of Section~\\ref{sec:inference}: a two-sided $p$-value and ",
                 "a 95\\%% interval obtained by inverting the randomization test over ",
                 "%s relabelings of which side of each boundary pair is rationed within ",
                 "each rotation cycle, and, on the line below, over the complete enumeration ",
                 "of all $2^{15}=32{,}768$ whole-boundary relabelings. Neither uses a normal ",
                 "approximation or any model of the spatial covariance. The control-side mean is"),
          format(length(By), big.mark = ",")), t1[nt], fixed = TRUE)
writeLines(t1, f1)

t4 <- c(
"% Generated by 12_randomization_inference.R -- do not edit by hand.",
"\\begin{table}[htbp]", "\\centering",
"\\caption{Exact design-based inference and diagnostic checks}",
"\\label{tab:designchecks}",
"\\begin{threeparttable}\\small",
"\\begin{tabular}{lcccc}", "\\toprule",
"Check & Statistic & Estimate & Null s.d. / SE & $p$-value \\\\",
"\\midrule",
"\\addlinespace[0.1cm]",
"\\multicolumn{5}{l}{\\textit{Panel A: Exact randomization inference on the pooled discontinuity}} \\\\",
"\\addlinespace[0.1cm]",
sprintf("Cycle-level relabeling (primary) & $\\hat\\tau_{\\mathrm{IV}}$ & $%s$ & %s & %s \\\\",
        fmt(obs_y), fmt(sd(By), 4), fmt(pB, 3)),
sprintf("\\quad 95\\%% randomization interval & & \\multicolumn{2}{c}{[$%s$, $%s$]} & \\\\",
        fmt(ciB[1]), fmt(ciB[2])),
sprintf("Whole-boundary relabeling (all $2^{15}$) & $\\hat\\tau_{\\mathrm{IV}}$ & $%s$ & %s & %s \\\\",
        fmt(obs_y), fmt(sd(Ay), 4), fmt(pA, 3)),
sprintf("\\quad 95\\%% randomization interval & & \\multicolumn{2}{c}{[$%s$, $%s$]} & \\\\",
        fmt(ciA[1]), fmt(ciA[2])),
sprintf("\\quad \\emph{memo}: meta-analytic SE & & & %s & \\\\", fmt(1/sqrt(WSUM), 4)),
sprintf("\\quad \\emph{memo}: date block bootstrap SE & & & %s & \\\\", fmt(0.008130, 4)),
"\\addlinespace[0.2cm]", "\\midrule", "\\addlinespace[0.1cm]",
"\\multicolumn{5}{l}{\\textit{Panel B: Heterogeneity}} \\\\",
"\\addlinespace[0.1cm]",
sprintf("Cochran's $Q$ ($\\chi^2_{%d}$ reference) & $Q$ & $%s$ & -- & %s \\\\",
        dfQ, formatC(Qstat, format="f", digits=2), fmt(pQ, 3)),
sprintf("Cochran's $Q$ (design-based reference) & $Q$ & $%s$ & -- & %s \\\\",
        formatC(Qstat, format="f", digits=2), fmt(pQ_ri, 3)),
sprintf("Familywise exact test, 15 boundaries & $\\max_b|z_b|$ & $%s$ & -- & %s \\\\",
        formatC(max(abs(Zobs)), format="f", digits=2), fmt(p_maxT, 3)),
"\\addlinespace[0.2cm]", "\\midrule", "\\addlinespace[0.1cm]",
"\\multicolumn{5}{l}{\\textit{Panel C: Further diagnostic checks}} \\\\",
"\\addlinespace[0.1cm]",
"Density test at the cutoff & \\citeauthor{cattaneo2018manipulation} $T$ & -- & -- & 0.066 \\\\",
"Sexual-violence comparison outcome & pooled $\\hat\\tau$ & $0.063$ & 0.088 & 0.474 \\\\",
"Treated-shift symmetry (4 of 9) & range $\\hat\\tau_s$ & \\multicolumn{2}{c}{$[-1.61,~0.63]$} & -- \\\\",
"Leave-one-boundary-out & range $\\hat\\tau_{\\mathrm{IV}}$ & \\multicolumn{2}{c}{$[-0.062,~-0.005]$} & -- \\\\",
"\\bottomrule", "\\end{tabular}",
"\\begin{tablenotes}[flushleft]\\small",
paste0("\\item \\textit{Notes:} Panel A reports exact design-based inference on the ",
"pooled inverse-variance weighted discontinuity, computed on the same balanced ",
"84,000-cell panel as Table~\\ref{tab:main_results}. The randomization exchanges which ",
"member of an adjacent polygon pair is the rationed one, which maps $\\hat\\tau_b$ to ",
"$-\\hat\\tau_b$ exactly and leaves the incidents, distances, bins, active dates and ",
"treated-day counts unchanged. The primary row randomises this within each nine-day ",
sprintf("rotation cycle (%s uniform draws from the relabeling group); the second row ",
        format(length(By), big.mark=",")),
"enumerates all $2^{15}=32{,}768$ whole-boundary relabelings and is the more conservative ",
"of the two. Intervals are obtained by inverting the randomization test, exploiting the ",
"linearity of the estimator in the cell counts; they involve no normal approximation and ",
"no model of the spatial covariance. Panel B tests $Q$ against both the $\\chi^2$ ",
"reference and the design-based one, and reports a Westfall--Young familywise exact test ",
"on the largest boundary-specific $z$-statistic. Panel C rows are computed on the earlier ",
"occupied-cell panel and are qualitative; see Section~\\ref{sec:designchecks}."),
"\\end{tablenotes}", "\\end{threeparttable}", "\\end{table}")
writeLines(t4, file.path(TAB, "table4.tex"))
cat(sprintf("\nwrote %s and patched %s\n", file.path(TAB,"table4.tex"), f1))
