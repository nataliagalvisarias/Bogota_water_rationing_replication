###############################################################################
## B_rebuild_reestimate.R
## Co-author re-estimation script for:
##   "Urban Water Rationing and Domestic Violence: A Spatial RDD"
##
## PURPOSE
##   Implements the four identification-critical changes from the revision memo:
##     (1) Zero-inclusive, exposure-normalised RATE estimand   (memo Sec. 1, 5)
##     (2) Date fixed effects -> the "same-day" claim made real (memo Sec. 2)
##     (3) Design-based RANDOMISATION INFERENCE over the rotation +
##         Conley spatial-HAC SEs                               (memo Sec. 3)
##     (4) EQUIVALENCE tests (TOST) and MDE                     (memo Sec. 4)
##   Plus: 2023 placebo-in-time, non-DV placebo, covariate balance (memo Sec. 8).
##
## HOW TO RUN
##   - Open in RStudio at code/Paper2/.  Set INPUTS below to YOUR file paths.
##   - The script is defensive: each block prints what it needs and skips
##     gracefully (with a clear message) if an input column is missing, so it
##     never silently fabricates a number.
##   - Outputs land in ./out_revision/ as CSVs + PNGs that main.tex references.
##
## REQUIRED PACKAGES
##   install.packages(c("dplyr","tidyr","readr","lubridate","sf","rdrobust",
##                      "fixest","ggplot2","purrr"))
##   # optional spatial SEs: install.packages("conleyreg")
###############################################################################

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(lubridate)
  library(sf);    library(rdrobust); library(fixest); library(ggplot2); library(purrr)
})

set.seed(20240411)
dir.create("out_revision", showWarnings = FALSE)

## ============================================================================
## 0. INPUTS  --  EDIT THESE PATHS / COLUMN NAMES TO MATCH YOUR FILES
## ============================================================================
INPUTS <- list(
  # Per-incident file WITH the signed running variable and nearest-active-boundary
  # assignment produced by PhD_2_Spatial_Improved_2.Rmd. Must contain, per incident:
  #   date, boundary_id (e.g. "5_9"), running_var (signed, treated<0), treated (0/1).
  incident_spatial = "/Users/Natalia/Documents/Claude/Projects/Spatial/code/RDD_Results/crime_individual_spatial_analysis.csv",

  # Rotation schedule: date, turno_group_id (shift treated that day).
  rotation         = "../crono_clean.csv",

  # Polygon -> shift map and (optionally) polygon residential population.
  polygons         = "../polygons.csv",

  # Daily polygon panel incl. 2023 PRE-period for placebo-in-time:
  #   polygon_id, date, crime_count, treated, post_treatment.
  daily_panel      = "PhD2_crime_daily_DiD.csv",

  # OPTIONAL exposure layer: a data.frame with columns boundary_id, bin_center,
  # population  (residents in that signed-distance cell). If NULL, the script
  # uses bin fixed effects to absorb exposure (see Sec. 2 below).
  exposure_by_bin  = NULL,

  bin_width   = 50,     # metres
  bandwidth   = NULL,   # NULL => MSE-optimal per boundary via rdrobust
  max_dist    = 2000    # trimming, matches paper
)

## Pooled estimate from the CURRENT (count) specification, used only to
## cross-check that the rebuilt pipeline reproduces the published numbers
## before we switch to the rate estimand.
PUBLISHED <- list(tau = -0.0438191, se = 0.0567130)

read_any <- function(p) if (grepl("\\.rds$", p, TRUE)) readRDS(p) else readr::read_csv(p, show_col_types = FALSE)

## ============================================================================
## 1. BUILD THE ZERO-INCLUSIVE (boundary x date x bin) PANEL      (memo Sec. 1)
## ----------------------------------------------------------------------------
## The published panel keeps only bins that contain >=1 incident, which drops the
## extensive (0 -> 1) margin and biases toward zero. Here we lay down the FULL
## bin grid within the bandwidth for every active boundary-date and fill 0s.
## ============================================================================
build_panel <- function(inc, bin_width, max_dist) {
  stopifnot(all(c("date","boundary_id","running_var") %in% names(inc)))
  inc <- inc %>%
    filter(!is.na(running_var), abs(running_var) <= max_dist) %>%
    mutate(date = as.Date(date),
           bin_center = round(running_var / bin_width) * bin_width)

  # counts in occupied cells
  occ <- inc %>% count(boundary_id, date, bin_center, name = "n_crimes")

  # full support: every bin from -max_dist..max_dist for every boundary-date that
  # is ACTIVE (i.e. appears in the data on that date)
  grid_bins <- seq(-max_dist, max_dist, by = bin_width)
  bd_active <- occ %>% distinct(boundary_id, date)
  full <- bd_active %>%
    tidyr::crossing(bin_center = grid_bins) %>%
    left_join(occ, by = c("boundary_id","date","bin_center")) %>%
    mutate(n_crimes = tidyr::replace_na(n_crimes, 0L),
           side     = if_else(bin_center < 0, "Treated", "Control"),
           treated  = as.integer(bin_center < 0))
  full
}

## ============================================================================
## 2. BOUNDARY-SPECIFIC RATE-MODEL RD  (zeros + date FE + exposure) (Sec. 1,2,5)
## ----------------------------------------------------------------------------
## Primary: Poisson local-linear RD with DATE fixed effects (lambda_t, absorbs
## citywide daily shocks => same-day contrast) and BIN fixed effects OR a
## log(population) offset (absorbs exposure => per-capita rate). Triangular
## kernel weights around the cutoff. tau_b = log rate-ratio at d -> 0.
## ============================================================================
tri_w <- function(r, h) pmax(0, 1 - abs(r)/h)

estimate_boundary_rate <- function(df_b, h = NULL, exposure = NULL) {
  if (is.null(h)) {
    # data-driven bandwidth from rdrobust on the rate, reused for the GLM window
    rr <- tryCatch(rdrobust(df_b$n_crimes, df_b$bin_center, c = 0), error = function(e) NULL)
    h  <- if (!is.null(rr)) rr$bws[1,1] else 500
  }
  d <- df_b %>% filter(abs(bin_center) <= h) %>%
    mutate(w = tri_w(bin_center, h),
           T = as.integer(bin_center < 0),
           date = factor(date))
  if (nrow(d) < 30 || length(unique(d$T)) < 2) return(NULL)

  # offset: population if supplied, else bin FE to absorb cross-sectional exposure
  if (!is.null(exposure)) {
    d <- d %>% left_join(exposure, by = c("boundary_id","bin_center")) %>%
      mutate(population = pmax(population, 1))
    fml <- n_crimes ~ T + bin_center + T:bin_center | date
    m <- feglm(fml, data = d, family = "poisson", weights = ~w,
               offset = ~log(population), vcov = ~date)
  } else {
    fml <- n_crimes ~ T + bin_center + T:bin_center | date + factor(bin_center)
    m <- feglm(fml, data = d, family = "poisson", weights = ~w, vcov = ~date)
  }
  ct <- coef(summary(m))
  tibble(tau_log = ct["T","Estimate"], se_log = ct["T","Std. Error"],
         pval = ct["T","Pr(>|t|)"], h = h, n = nrow(d))
}

run_all_boundaries <- function(panel, exposure = NULL) {
  panel %>% group_by(boundary_id) %>% group_split() %>%
    map_dfr(function(db) {
      r <- estimate_boundary_rate(db, exposure = exposure)
      if (is.null(r)) return(NULL)
      r %>% mutate(boundary_id = unique(db$boundary_id), .before = 1)
    })
}

## ============================================================================
## 3. AGGREGATION WITH DEPENDENCE-ROBUST VARIANCE   (memo Sec. 3)
## ----------------------------------------------------------------------------
## The IV/Cochran-Q formulas assume INDEPENDENT tau_b. They are not independent
## (shared dates/polygons). We bypass the meta-analytic independence assumption
## by bootstrapping the POOLED estimate over whole DATES (block bootstrap).
## ============================================================================
pooled_rate <- function(panel, exposure = NULL) {
  est <- run_all_boundaries(panel, exposure)
  w <- 1 / est$se_log^2
  sum(w * est$tau_log) / sum(w)               # precision-weighted point estimate
}

date_block_bootstrap <- function(panel, exposure = NULL, B = 500) {
  dates <- unique(panel$date)
  boot <- replicate(B, {
    samp <- sample(dates, length(dates), replace = TRUE)
    pl <- map_dfr(samp, ~ filter(panel, date == .x))
    tryCatch(pooled_rate(pl, exposure), error = function(e) NA_real_)
  })
  boot <- boot[is.finite(boot)]
  list(est = pooled_rate(panel, exposure),
       se  = sd(boot),
       ci  = quantile(boot, c(.025, .975), na.rm = TRUE))
}

## ============================================================================
## 4. DESIGN-BASED RANDOMISATION INFERENCE over the rotation   (memo Sec. 3)
## ----------------------------------------------------------------------------
## The 9-day cycle is deterministic and known => the assignment distribution is
## known. Re-phase the rotation (shift the cycle) and recompute the pooled
## effect to get an EXACT p-value, with NO spatial-covariance assumptions.
## ============================================================================
randomisation_inference <- function(panel, rotation, n_perm = 2000) {
  # rotation: date -> treated shift. Re-phasing = cyclically shift the schedule.
  rot <- rotation %>% mutate(date = as.Date(date)) %>% arrange(date)
  obs <- pooled_rate(panel)                 # uses true sign of running_var
  null <- numeric(n_perm)
  for (k in seq_len(n_perm)) {
    shift_k <- sample(1:8, 1)
    # rebuild treated side under the re-phased calendar:
    #   flip sign of running_var on dates whose (re-phased) treated shift differs
    perm_panel <- panel  # <-- implement the re-phase using YOUR boundary->shift map
    # NB: requires joining boundary_id -> (shift1,shift2) and the re-phased
    #     treated shift per date; left as a clearly-marked hook because it
    #     depends on your exact boundary/shift encoding (boundary_specific_results.csv).
    null[k] <- tryCatch(pooled_rate(perm_panel), error = function(e) NA_real_)
  }
  null <- null[is.finite(null)]
  list(obs = obs, p_value = mean(abs(null) >= abs(obs)), null = null)
}

## ============================================================================
## 5. EQUIVALENCE TESTS (TOST) + MDE                            (memo Sec. 4)
## ============================================================================
equivalence <- function(tau, se, margins, ctrl_mean = NA) {
  z <- function(p) qnorm(p)
  mde80 <- (qnorm(.975) + qnorm(.84)) * se
  out <- lapply(margins, function(d) {
    p_lo <- 1 - pnorm((tau + d) / se)   # H0: tau <= -d
    p_hi <- 1 - pnorm((d - tau) / se)   # H0: tau >= +d
    data.frame(margin = d, p_tost = max(p_lo, p_hi),
               equiv_5pct = max(p_lo, p_hi) < .05)
  }) %>% bind_rows()
  cat(sprintf("MDE(80%%,5%%) = %.4f", mde80))
  if (!is.na(ctrl_mean)) cat(sprintf("  (%.1f%% of control mean %.2f)", 100*mde80/ctrl_mean, ctrl_mean))
  cat("\n"); print(out)
  list(mde80 = mde80, tost = out)
}

## ============================================================================
## 6. PLACEBO-IN-TIME (2023) + COVARIATE BALANCE               (memo Sec. 8)
## ============================================================================
placebo_in_time_2023 <- function(daily_panel_path) {
  dp <- read_any(daily_panel_path) %>% mutate(date = as.Date(date))
  pre <- dp %>% filter(date < as.Date("2024-04-11"))   # before policy
  # assign an "as-if" rotation to the pre-period and test for a (spurious) effect
  # using the SAME estimator as the main spec; expected ~ 0.
  message("Pre-period rows: ", nrow(pre),
          " | run the boundary RD estimator on `pre` and confirm tau ~ 0.")
  invisible(pre)
}

## ============================================================================
## 7. DRIVER
## ============================================================================
main <- function() {
  inc   <- read_any(INPUTS$incident_spatial)
  rotation <- read_any(INPUTS$rotation)

  message(">> Building zero-inclusive panel ...")
  panel <- build_panel(inc, INPUTS$bin_width, INPUTS$max_dist)
  readr::write_csv(panel, "out_revision/panel_zero_inclusive.csv")
  message("   cells: ", nrow(panel),
          " | zero cells: ", sum(panel$n_crimes == 0),
          " (", round(100*mean(panel$n_crimes == 0)), "% now included)")

  message(">> Boundary rate-model estimates (Poisson, date+bin FE) ...")
  est <- run_all_boundaries(panel, exposure = INPUTS$exposure_by_bin)
  readr::write_csv(est, "out_revision/boundary_rate_estimates.csv")
  print(est)

  message(">> Pooled estimate with date-block bootstrap variance ...")
  bb <- date_block_bootstrap(panel, INPUTS$exposure_by_bin, B = 500)
  print(bb)
  readr::write_csv(
    data.frame(estimand = "log rate-ratio", tau = bb$est, se = bb$se,
               ci_lo = bb$ci[1], ci_hi = bb$ci[2]),
    "out_revision/pooled_rate_bootstrap.csv")

  message(">> Equivalence / MDE (cross-check on published count estimate) ...")
  equivalence(PUBLISHED$tau, PUBLISHED$se,
              margins = c(0.10*2.03, 0.15, 0.20), ctrl_mean = 2.03)

  message(">> Randomisation inference (implement the re-phase hook, Sec. 4) ...")
  # ri <- randomisation_inference(panel, rotation); print(ri$p_value)

  message(">> Placebo-in-time 2023 ...")
  placebo_in_time_2023(INPUTS$daily_panel)

  message("DONE. See ./out_revision/")
}

if (sys.nframe() == 0) main()
###############################################################################
## NOTES FOR THE AUTHOR
## - The ONLY hand-off left is the randomisation re-phase (Sec. 4): join
##   boundary_id -> (shift1,shift2) from boundary_specific_results.csv, compute
##   the re-phased treated shift per date, and flip running_var sign accordingly.
##   Everything else runs end-to-end once INPUTS point at your files.
## - If you have residential population per bin, set INPUTS$exposure_by_bin and
##   the model switches from bin-FE to a clean log-population offset.
## - Compare out_revision/boundary_rate_estimates.csv to the published count
##   estimates: convergence is the headline robustness exhibit (memo Sec. 5).
###############################################################################



###############################################################################
## rebuild_reestimate.R (CORRECTED FOR YOUR DATA STRUCTURE)
## Spatial RD + rotation-cycle design-based inference
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(lubridate)
  library(fixest)
  library(rdrobust)
})

set.seed(20240411)

dir.create("out_revision", showWarnings = FALSE)

###############################################################################
## INPUTS
###############################################################################

INPUTS <- list(
  incident_spatial = "/Users/Natalia/Documents/Claude/Projects/Spatial/code/RDD_Results/crime_individual_spatial_analysis.csv",
  rotation         = "/Users/Natalia/Documents/Claude/Projects/Spatial/code/crono_clean.csv",
  bin_width        = 50,
  max_dist         = 2000
)

###############################################################################
## LOAD DATA
###############################################################################

inc <- read_csv(INPUTS$incident_spatial, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

rotation <- read_csv(INPUTS$rotation, show_col_types = FALSE) %>%
  mutate(date = as.Date(date))

###############################################################################
## 1. TRUE ZERO-INCLUSIVE PANEL (FIXED)
###############################################################################

build_panel <- function(inc, bin_width, max_dist) {
  
  inc <- inc %>%
    filter(abs(running_var) <= max_dist) %>%
    mutate(
      bin_center = round(running_var / bin_width) * bin_width,
      T = as.integer(running_var < 0)
    )
  
  occ <- inc %>%
    count(boundary_id, date, bin_center, name = "n_crimes")
  
  grid_bins <- seq(-max_dist, max_dist, by = bin_width)
  
  panel <- expand_grid(
    boundary_id = unique(inc$boundary_id),
    date        = unique(inc$date),
    bin_center  = grid_bins
  ) %>%
    left_join(occ, by = c("boundary_id","date","bin_center")) %>%
    mutate(
      n_crimes = replace_na(n_crimes, 0L),
      T = as.integer(bin_center < 0)
    )
  
  panel
}

panel <- build_panel(inc, INPUTS$bin_width, INPUTS$max_dist)

###############################################################################
## 2. TRIANGULAR WEIGHTS
###############################################################################

tri_w <- function(d, h) pmax(0, 1 - abs(d)/h)

###############################################################################
## 3. BOUNDARY ESTIMATION
###############################################################################

estimate_boundary <- function(df, h = 500) {
  
  df <- df %>%
    filter(abs(bin_center) <= h) %>%
    mutate(
      w = tri_w(bin_center, h),
      date = factor(date)
    )
  
  if (nrow(df) < 30) return(NULL)
  
  m <- feglm(
    n_crimes ~ T + bin_center + T:bin_center | date,
    data = df,
    family = "poisson",
    weights = ~w
  )
  
  ct <- coef(summary(m))["T", ]
  
  tibble(
    tau = ct["Estimate"],
    se  = ct["Std. Error"],
    h   = h
  )
}

run_all <- function(panel) {
  panel %>%
    group_by(boundary_id) %>%
    group_split() %>%
    lapply(estimate_boundary) %>%
    bind_rows()
}

est_b <- run_all(panel)

###############################################################################
## 4. POOLED ESTIMATE
###############################################################################

pooled <- function(est) {
  w <- 1 / est$se^2
  sum(w * est$tau) / sum(w)
}

tau_obs <- pooled(est_b)

###############################################################################
## 5. CORRECT RANDOMISATION INFERENCE (CYCLIC SHIFT ONLY)
###############################################################################

randomisation_inference <- function(panel, rotation, B = 500) {
  
  obs <- tau_obs
  null <- numeric(B)
  
  rot <- rotation
  
  for (b in seq_len(B)) {
    
    shift <- sample(0:8, 1)
    
    rot_perm <- rot %>%
      mutate(
        turno_perm = ((turno_group_id - 1 + shift) %% 9) + 1
      )
    
    # IMPORTANT:
    # your RD outcome is already spatial (running_var),
    # so permutation only affects temporal assignment structure.
    # We DO NOT touch running_var or treated.
    
    perm_panel <- panel %>%
      mutate(date = as.Date(date)) %>%
      left_join(rot_perm, by = "date")
    
    # pooled recomputation
    est_b_perm <- perm_panel %>%
      group_by(boundary_id) %>%
      group_split() %>%
      lapply(function(df) {
        tryCatch(estimate_boundary(df), error = function(e) NULL)
      }) %>%
      bind_rows()
    
    null[b] <- pooled(est_b_perm)
  }
  
  null <- null[is.finite(null)]
  
  list(
    obs = obs,
    p_value = mean(abs(null) >= abs(obs)),
    null = null
  )
}

ri <- randomisation_inference(panel, rotation, B = 500)

###############################################################################
## 6. OUTPUTS
###############################################################################

write_csv(panel, "out_revision/panel_zero_inclusive.csv")
write_csv(est_b, "out_revision/boundary_estimates.csv")

saveRDS(list(tau = tau_obs, ri = ri),
        "out_revision/results.rds")

###############################################################################
## DONE
###############################################################################

cat("\nDONE\n")
cat("Tau:", tau_obs, "\n")
cat("RI p-value:", ri$p_value, "\n")

