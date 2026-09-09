#!/usr/bin/env Rscript
## =============================================================================
## 11_final_figures_and_checks.R
## Figures 3-7 and the bandwidth / placebo checks, all computed on the SAME
## balanced zero-filled panel that produces Tables 1-3, so that no exhibit in
## the paper rests on a different sample from the one it illustrates.
## Run after 10_final_submission_estimates.R.
## =============================================================================
suppressPackageStartupMessages({library(dplyr); library(readr); library(ggplot2); library(tidyr)})
ROOT <- Sys.getenv("ROOT", "."); DER <- file.path(ROOT,"data","derived")
FIG  <- file.path(ROOT,"figures"); OUT <- file.path(ROOT,"output")
dir.create(FIG, recursive=TRUE, showWarnings=FALSE)
PURPLE <- "#660099"; LILAC <- "#B280CC"
panel <- read_csv(file.path(DER,"panel_balanced_zerofilled.csv"), show_col_types=FALSE)
est   <- read_csv(file.path(OUT,"boundary_estimates_final.csv"),  show_col_types=FALSE)
bws   <- read_csv(file.path(DER,"bandwidths_mserd.csv"),          show_col_types=FALSE)
S     <- read_csv(file.path(OUT,"manuscript_statistics.csv"),     show_col_types=FALSE)
g <- function(k) S$value[S$quantity==k]
BASELINE <- g("control_baseline_per_cell"); TAU <- g("pooled_tau"); SE <- g("pooled_se_bootstrap")

cluster_vcov <- function(X,u,w,cl){
  bread <- solve(crossprod(X*sqrt(w))); meat <- matrix(0,ncol(X),ncol(X))
  for (gg in split(seq_along(u), cl)) { s <- crossprod(X[gg,,drop=FALSE], w[gg]*u[gg]); meat <- meat + tcrossprod(s) }
  G <- length(unique(cl)); N <- length(u); K <- ncol(X)
  bread %*% meat %*% bread * (G/(G-1)) * ((N-1)/(N-K))
}
rd_ll <- function(y,x,h,cl){
  w <- pmax(0,1-abs(x)/h); k <- w>0
  if (sum(k)<8L) return(NULL)
  y<-y[k];x<-x[k];w<-w[k];cl<-cl[k]; D<-as.numeric(x<0)
  X <- cbind(1,D,x,x*D); if (qr(X*sqrt(w))$rank<ncol(X)) return(NULL)
  b <- solve(crossprod(X*sqrt(w)), crossprod(X,w*y)); u <- as.vector(y-X%*%b)
  V <- cluster_vcov(X,u,w,cl); list(tau=unname(b[2]), se=sqrt(V[2,2]))
}
ivw <- function(t,s){ w<-1/s^2; tt<-sum(w*t)/sum(w); list(tau=tt, se=1/sqrt(sum(w))) }
pool_at <- function(p, hfun){
  r <- lapply(sort(unique(p$boundary_id)), function(b){
    gg <- p[p$boundary_id==b,]; f <- rd_ll(gg$n, gg$bin, hfun(b), as.character(gg$date))
    if (is.null(f)) NULL else data.frame(tau=f$tau, se=f$se)}) |> bind_rows()
  if (!nrow(r)) return(NULL); ivw(r$tau, r$se)
}

## ---- Figure 3: forest plot ---------------------------------------------------
fp <- est |> mutate(lab=gsub("_","--",boundary_id)) |> arrange(tau) |> mutate(lab=factor(lab,levels=lab))
p3 <- ggplot(fp, aes(tau, lab)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey45") +
  geom_vline(xintercept=TAU, colour=PURPLE, linewidth=0.9) +
  geom_errorbarh(aes(xmin=ci_lower, xmax=ci_upper), height=0, colour=PURPLE, alpha=.75) +
  geom_point(shape=21, fill="white", colour=PURPLE, size=2.4, stroke=.9) +
  labs(x=expression(hat(tau)[b]~"(incidents per boundary "%*%" date "%*%" bin cell)"),
       y="Boundary (polygon pair)") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"Figure_3.pdf"), p3, width=7, height=5, device=cairo_pdf)

## ---- Figure 4: histogram -----------------------------------------------------
p4 <- ggplot(est, aes(tau)) +
  geom_histogram(binwidth=.05, fill=LILAC, colour="white", boundary=0) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey45") +
  geom_vline(xintercept=TAU, colour=PURPLE, linewidth=1) +
  labs(x=expression("Boundary-specific estimate"~hat(tau)[b]), y="Number of boundaries") +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"Figure_4.pdf"), p4, width=6, height=4.2, device=cairo_pdf)

## ---- Figure 5: power curve and equivalence bounds ----------------------------
eff <- seq(-0.05, 0.05, length.out=400)
pw  <- pnorm(abs(eff)/SE - qnorm(0.975)) + pnorm(-abs(eff)/SE - qnorm(0.975))
mde <- 2.80*SE; hi <- TAU + 1.96*SE
p5 <- ggplot(data.frame(eff, pw), aes(eff, pw)) +
  geom_hline(yintercept=.8, linetype="dotted", colour="grey45") +
  geom_line(colour=PURPLE, linewidth=.9) +
  geom_vline(xintercept=mde, linetype="dashed", colour=PURPLE) +
  geom_vline(xintercept=hi, linetype="dotdash", colour="grey35") +
  scale_x_continuous(name="True effect (incidents per cell)",
    sec.axis=sec_axis(~ 100*./BASELINE, name="Percent of the control-side baseline")) +
  scale_y_continuous(name="Power", limits=c(0,1)) +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"Figure_5.pdf"), p5, width=6.5, height=4.2, device=cairo_pdf)

## ---- Figure 6: bandwidth sensitivity -----------------------------------------
hs <- seq(150, 1500, by=75)
bw_sens <- lapply(hs, function(h){
  r <- pool_at(panel, function(b) h); if (is.null(r)) return(NULL)
  data.frame(h=h, tau=r$tau, lo=r$tau-1.96*r$se, hi=r$tau+1.96*r$se)}) |> bind_rows()
write_csv(bw_sens, file.path(OUT,"bandwidth_sensitivity_final.csv"))
p6 <- ggplot(bw_sens, aes(h, tau)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey45") +
  geom_vline(xintercept=mean(est$bandwidth), linetype="dashed", colour=PURPLE, alpha=.7) +
  geom_linerange(aes(ymin=lo, ymax=hi), colour=PURPLE, alpha=.8) +
  geom_point(colour=PURPLE, size=1.9) +
  labs(x="Bandwidth (metres from the boundary)",
       y=expression("Pooled"~hat(tau)[IV]~"(incidents per cell)")) +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"Figure_6.pdf"), p6, width=6.5, height=4.2, device=cairo_pdf)

## ---- Figure 7: placebo cutoffs ------------------------------------------------
## Placebo cutoffs must not see the true boundary. A window centred on a
## displaced cutoff that straddles zero picks up the real discontinuity, which
## is why naive displaced-cutoff placebos over-reject. Following standard
## practice, each placebo uses ONLY the observations on one side of the true
## cutoff -- the side the artificial cutoff sits on -- and a window symmetric
## about the artificial cutoff, with the bandwidth capped by its half-width.
shifts <- c(-750,-500,-250,250,500,750)
plac_one <- lapply(shifts, function(d){
  half <- min(abs(d), 2000 - abs(d))
  p <- panel |> filter(sign(bin) == sign(d)) |>
       mutate(bin = bin - d) |> filter(abs(bin) <= half)
  if (!nrow(p)) return(NULL)
  r <- pool_at(p, function(b) min(bws$h_mserd[match(b,bws$boundary_id)], half))
  if (is.null(r)) return(NULL)
  data.frame(shift_m=d, half_window_m=half, tau=r$tau, se=r$se,
             lo=r$tau-1.96*r$se, hi=r$tau+1.96*r$se,
             pval=2*pnorm(-abs(r$tau/r$se)))}) |> bind_rows()
plac <- bind_rows(
  data.frame(shift_m=0, half_window_m=2000, tau=TAU, se=SE,
             lo=TAU-1.96*SE, hi=TAU+1.96*SE, pval=2*pnorm(-abs(TAU/SE))),
  plac_one) |> arrange(shift_m)
write_csv(plac, file.path(OUT,"placebo_cutoffs_final.csv"))
p7 <- ggplot(plac, aes(shift_m, tau)) +
  geom_hline(yintercept=0, linetype="dotted", colour="grey45") +
  geom_linerange(aes(ymin=lo, ymax=hi), colour=PURPLE, alpha=.8) +
  geom_point(aes(shape=shift_m==0), colour=PURPLE, size=2.6, fill="white", stroke=.9) +
  scale_shape_manual(values=c(`TRUE`=19,`FALSE`=21), guide="none") +
  labs(x="Displacement of the cutoff (metres)",
       y=expression("Pooled"~hat(tau)[IV]~"(incidents per cell)")) +
  theme_minimal(base_size=11) + theme(panel.grid.minor=element_blank())
ggsave(file.path(FIG,"Figure_7.pdf"), p7, width=6.5, height=4.2, device=cairo_pdf)

cat("\n--- bandwidth sensitivity (pooled) ---\n"); print(as.data.frame(bw_sens), digits=3, row.names=FALSE)
cat("\n--- placebo cutoffs (pooled) ---\n");        print(as.data.frame(plac),    digits=3, row.names=FALSE)
cat("\nFigures 3-7 written to", FIG, "\n")

## ---- Placebo-calibrated conservative inference --------------------------------
## The six one-sided placebo estimates are draws from a distribution that ought
## to be centred on zero. Their spread is an empirical benchmark for how much
## the pooled estimator moves when nothing is happening, and it is wider than
## the nominal standard error. Reporting equivalence bounds against BOTH the
## nominal and the placebo-calibrated scale is the conservative course.
pl <- plac[plac$shift_m != 0, ]
SE_PLACEBO <- sd(pl$tau)
calib <- data.frame(
  basis = c("date-block bootstrap", "placebo-calibrated (conservative)"),
  se    = c(SE, SE_PLACEBO),
  mde_pct       = 100 * 2.80 * c(SE, SE_PLACEBO) / BASELINE,
  upper_ci_pct  = 100 * (TAU + 1.96 * c(SE, SE_PLACEBO)) / BASELINE,
  tightest_pct  = 100 * (abs(TAU) + qnorm(0.95) * c(SE, SE_PLACEBO)) / BASELINE,
  tost_p        = pmax(pnorm((TAU - 0.10*BASELINE)/c(SE,SE_PLACEBO)),
                       pnorm((TAU + 0.10*BASELINE)/c(SE,SE_PLACEBO), lower.tail = FALSE)))
write_csv(calib, file.path(OUT, "equivalence_calibrated.csv"))
cat("\n--- equivalence bounds, nominal vs placebo-calibrated ---\n")
print(calib, digits = 3, row.names = FALSE)
cat(sprintf("\nplacebo spread SD = %.5f  (%.1fx the bootstrap SE of %.5f)\n",
            SE_PLACEBO, SE_PLACEBO/SE, SE))
cat(sprintf("|tau| at the true cutoff = %.5f; rank among the 7 cutoffs by |tau| = %d of 7\n",
            abs(TAU), rank(abs(plac$tau))[plac$shift_m == 0]))
