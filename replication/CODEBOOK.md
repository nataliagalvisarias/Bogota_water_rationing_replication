# Codebook

Replication package for **Does Urban Water Scarcity Increase Domestic Violence? Evidence
from Water Rationing in Bogotá** (Journal of Urban Economics).

All spatial data are in **EPSG:3116** (MAGNA-SIRGAS / Colombia Bogotá zone); distances and
coordinates are in **metres**. Dates are `YYYY-MM-DD`.

---

## Conventions that govern every file

**Sign of the running variable and of the treatment effect.** Treated observations lie on
the **negative** side of the running variable (equation 6 of the manuscript). The
treatment effect is defined as the treated-side limit minus the control-side limit, so a
**positive `tau` means rationing raises reported violence**. This is the opposite of what
`rdrobust` returns: that package reports the right limit minus the left limit, which under
this sign convention is control minus treated. `10_final_submission_estimates.R` estimates
the treatment effect directly and does not need to be sign-flipped; if you call `rdrobust`
yourself on these files, negate its coefficient before interpreting it.

**Bins are half-open and never span the cutoff.** Treated side `[-50,0)`, `[-100,-50)`, …
with midpoints −25, −75, …; control side `[0,50)`, `[50,100)`, … with midpoints +25, +75, ….
Rounding the running variable to the nearest multiple of 50 would create a cell covering
`[-25,+25)` that mixes treated and control observations at the discontinuity itself.

**The estimation panel is balanced and zero-filled.** For every boundary and every date on
which that boundary separates a treated from an untreated polygon, all 80 bins spanning
±2,000 m are present, and cells with no incident carry `n = 0`. Counting only occupied
cells conditions the estimator on a positive outcome and attenuates the discontinuity
toward zero.

---

## `data/derived/rd_estimation_sample.csv` — incident level

20,766 rows before screening. Filter to `in_estimation_sample == 1 & is_duplicate_record == 0`
for the 20,527 incidents used in the paper.

| Variable | Type | Definition |
|---|---|---|
| `crime_id` | character | Incident identifier from the 123 Emergency Line. Not unique across rows; see `is_duplicate_record`. |
| `date` | date | Calendar date of the report. |
| `polygon_id` | integer | EAAB hydraulic polygon containing the incident. Identifiers run {1,2,3,5,…,13}; there is no polygon 4 in the rationing system. |
| `shift_id` | integer | Rationing shift (1–9) of that polygon. |
| `treated` | integer | 1 if the polygon was under a 24-hour suspension on that date. |
| `boundary_id` | character | Assigned treatment–control boundary, lower polygon id first, e.g. `3_11`. |
| `boundary_distance` | numeric | Unsigned distance to the assigned boundary, metres (0 to 2,000). |
| `running_var` | numeric | Signed running variable, metres. **Negative on the treated side.** |
| `treated_polygon_id`, `control_polygon_id` | integer | The two sides of the assigned boundary on that date. |
| `x_magna_m`, `y_magna_m` | numeric | Easting and northing, EPSG:3116, metres. **Not** degrees; the staging file called these `longitude`/`latitude`, which was wrong. |
| `is_duplicate_record` | integer | 1 if the row exactly repeats another on (`crime_id`, `date`, `boundary_id`). 202 such rows, an artefact of a many-to-many join to the adjacency table. Flagged rather than deleted so both the published and corrected sample sizes are recoverable. |
| `in_estimation_sample` | integer | 1 for the 15 estimation boundaries; 0 for boundary `1_8`, which has 37 incidents and cannot support a local fit. |
| `within_distance_screen` | integer | 1 if within ±2,000 m. |
| `abs_distance_m`, `log_distance_m`, `distance_km` | numeric | Convenience transforms of `running_var`. |

## `data/derived/panel_balanced_zerofilled.csv` — estimation panel

84,000 rows = 1,050 active boundary-dates × 80 bins. Written by step 10; the primary
estimation input. 83.1% of rows have `n = 0`.

| Variable | Type | Definition |
|---|---|---|
| `boundary_id` | character | Treatment–control boundary. |
| `date` | date | Calendar date on which the boundary is active. |
| `bin` | numeric | Bin midpoint in metres, one of ±25, ±75, …, ±1975. Negative = treated side. |
| `n` | integer | Incidents in the cell. Zero where none was reported; these rows are the point of the file. |

Control-side mean of `n` is **0.281**, the baseline against which every percentage in
Section 5 is expressed.

## `data/derived/bandwidths_mserd.csv` — MSE-optimal bandwidths

One row per boundary. Produced by `rdrobust` (Calonico–Cattaneo–Titiunik) on the balanced
panel and shipped so that step 10 is reproducible where `rdrobust` cannot be installed.
When the package *is* available, step 10 recomputes these and overwrites the file.

| Variable | Definition |
|---|---|
| `boundary_id` | Boundary. |
| `h_mserd` | MSE-optimal main bandwidth, metres (335–797, mean 504). |
| `b_mserd` | Bias bandwidth, metres. |
| `tau_rdrobust`, `se_rdrobust` | `rdrobust`'s own coefficient and robust SE, **on the package's control-minus-treated convention**. Retained for cross-checking only. |

## `output/boundary_estimates_final.csv` — Table 1, Panel A

| Variable | Definition |
|---|---|
| `boundary_id` | Boundary. |
| `tau`, `se`, `z`, `pval`, `ci_lower`, `ci_upper` | Treatment effect (treated − control) and date-clustered inference. |
| `bandwidth` | MSE-optimal bandwidth, metres. |
| `n_cells` | Cells for the boundary in the balanced panel. |
| `n_eff`, `n_eff_left`, `n_eff_right` | Cells with positive kernel weight, total and by side. |
| `n_incidents` | Incidents within 2,000 m of the boundary. |
| `n_dates` | Dates on which the boundary is active — 67 to 72, **not** the 325-day study window. |

## `output/manuscript_statistics.csv`

Every scalar quoted in the manuscript, as `quantity,value`. If a number appears in the
text it is in this file; nothing is typed by hand.

## Other outputs

| File | Contents |
|---|---|
| `output/equivalence_calibrated.csv` | Equivalence bounds on the bootstrap and the placebo-calibrated standard error. |
| `output/placebo_cutoffs_final.csv` | One-sided placebo estimates at ±250, ±500, ±750 m. |
| `output/bandwidth_sensitivity_final.csv` | Pooled estimate at bandwidths from 150 to 1,500 m. |
| `tables/table1.tex` … `table_c1.tex` | Production LaTeX, `\input` directly by the manuscript. |
| `figures/Figure_3.pdf` … `Figure_7.pdf` | Vector figures regenerated from the balanced panel. |

---

## Cleaning decisions, in the order applied

1. **Adjacency screen** — untreated polygons not contiguous to a treated polygon on that
   date are dropped (87,298 observations, 63.2% of the extract).
2. **Distance screen** — incidents beyond 2,000 m of the nearest active boundary dropped.
3. **De-duplication** — 202 exact repeats removed; N falls from 20,729 to **20,527**.
4. **Boundary screen** — `1_8` excluded for insufficient near-boundary support (37 incidents).
5. **Half-open binning** — no cell spans the cutoff.
6. **Grid completion** — all 80 bins per active boundary-date, zeros filled.

Steps 3, 5 and 6 are the ones that move the estimate; Section 5.3 of the manuscript
quantifies each. The full transformation log for the upstream curation is
`claude_workspace/curated_data/_curation_log.csv`.

## Software

R ≥ 4.1 with `dplyr`, `tidyr`, `readr`, `ggplot2`, `sf`; `rdrobust` optional (see above).
Step 12 additionally uses `ggplot2` and `cairo_pdf`; it allocates two 199,999 × 15 numeric
matrices, so allow roughly 100 MB of memory and about 80 seconds of runtime.
Python 3.10+ with `numpy`, `pandas`, `geopandas`, `matplotlib` for `06_design_checks.py`
and `08_python_reconstruction.py`.

## Reproducing

```bash
Rscript code/00_master.R      # runs steps 01, 04, 05, 07, 10, 11, 12
python3 code/06_design_checks.py   # density test and comparison outcomes only
```

Every script prints its recomputed values against those in the manuscript and flags any
divergence with `***`. A clean run prints none.


## Step 12 — exact design-based inference

`code/12_randomization_inference.R` is the source of the paper's primary inferential
result. It reads `data/derived/rd_estimation_sample.csv` and
`data/derived/panel_balanced_zerofilled.csv`, asserts that the two agree cell for cell, and
computes an exact randomization test of the pooled discontinuity.

**Randomization used.** Which member of an adjacent polygon pair is the rationed one is
relabeled. This maps τ̂_b to −τ̂_b exactly and preserves the incidents, their distances and
bins, the active dates, the panel dimensions and the number of rationed days on each side.
Cyclic re-phasing of the rationing calendar — the construction used in earlier drafts — is
**not** valid on this sample and is not used; the header of the script documents the
diagnostic in full.

| Level | Group | Draws | Role |
|---|---|---|---|
| boundary × nine-day cycle | 2^(Σ n_b^cyc) | 199,999 uniform | primary |
| whole boundary | 2^15 = 32,768 | complete enumeration | conservative |

**Closed form.** The estimator is linear in the cell counts, so
τ̂_b(f) = τ̂_b(0) + f′v_b with v_b = (A0_b − A1_b)(aN_b − aP_b), where a_b is the row of
(n_b G_b)^−1 U_b returning the discontinuity coefficient and A1_b, A0_b are the
block-aggregated count matrices on the realised treated and control sides. Confidence
intervals come from inverting the test: under H: τ = τ₀ the relabeled statistic is exactly
T_g(y) − τ₀ T_g(D).

**Assertions.** The script halts unless the reconstruction reproduces the 84,000-cell panel
exactly, the fifteen boundary estimates match `output/boundary_estimates_final.csv` to
1e-9, a full relabeling maps every τ̂_b to −τ̂_b to 1e-9, T(D) at the identity equals 1 to
1e-8, and each inverted interval is strictly interior to its search grid.

**Outputs.**

| File | Contents |
|---|---|
| `output/randomization_inference.csv` | 38 named quantities: *p*-values, intervals, equivalence bounds, heterogeneity statistics, policy magnitudes |
| `output/ri_boundary_pvalues.csv` | per-boundary τ̂, SE, *z*, exact *p*, incident count |
| `output/ri_null_distribution.csv` | the 32,767 enumerated null draws and their T(D) values |
| `tables/table4.tex` | Table 4, design-based inference and diagnostics |
| `tables/table1.tex` | Panels B and C patched with the design-based rows (idempotent) |
| `figures/Figure_8.pdf` | both null distributions |

**Environment variable.** `CITYWIDE_DV_PER_DAY` (default 302.8) sets the citywide daily
domestic-violence baseline used for the policy magnitudes. Nothing else in the paper
depends on it.
