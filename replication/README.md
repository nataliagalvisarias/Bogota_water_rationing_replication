# Replication package

**Does Urban Water Scarcity Increase Domestic Violence? Evidence from Water Rationing in Bogotá**

Natalia Galvis Arias (University of Manchester)
Prepared for the *Journal of Urban Economics*.

This README follows the [template for social science replication packages](https://social-science-data-editors.github.io/template_README/) used by the economics data editors.

---

## 1. Overview

The code in this replication package constructs the estimation sample and produces all tables and figures in the paper from three sources: the hydraulic service polygon layer published by the Bogotá Water and Sewerage Company (EAAB), the rationing calendar established by Resolución 291 de 2024, and incident-level emergency-call records held by the Bogotá Secretariat of Security.

**The incident-level records are restricted and are not included in this package.** Section 4 documents how to obtain them. To keep the results verifiable without them, the package ships the derived boundary-level estimates and the pooled distance-bin panel from which every reported number is computed. A replicator with no data access can therefore reproduce every table and figure in the paper; a replicator with data access can additionally rebuild those derived files from the raw records.

Two execution paths are supported:

| Path | Requires restricted data | Scripts run | What it reproduces | Runtime |
|---|---|---|---|---|
| **A. Derived-data (default)** | No | 01, 04, 05, 07 | Every table and figure in the paper | < 2 minutes |
| **B. Full** | Yes | 01–07 | The above, plus the derived files themselves | 20–40 minutes |

---

## 2. Data availability statement

### 2.1 Rights

The author certifies that she has legitimate access to, and permission to use, all data employed in the manuscript. The author certifies that she has permission to redistribute the public and derived data included here, and does **not** have permission to redistribute the restricted incident-level records.

### 2.2 Summary of data sources

| Dataset | Source | Provided here | Licence |
|---|---|---|---|
| 123 Emergency Line incident records (crimes against women, Jan 2023 – May 2025) | Secretaría Distrital de Seguridad, Convivencia y Justicia, Bogotá | **No** — restricted | Data-access agreement; redistribution not permitted |
| Hydraulic service polygons (12 polygons, EPSG:3116) | EAAB, via [mapas.bogota.gov.co](https://mapas.bogota.gov.co/?l=72004) | Yes — `data/public/hydraulic_polygons.gpkg` | Open government data |
| Rationing calendar, 11 Apr 2024 – 1 Apr 2025 (325 days) | Constructed from Resolución 291 de 2024 and Resolución 980 de 2024 | Yes — `data/public/rationing_calendar.csv` | Public administrative act |
| Polygon-to-shift crosswalk | Resolución 291 de 2024 (Appendix Table A.1) | Yes — `data/public/polygon_shift_crosswalk.csv` | Public administrative act |
| Boundary-level RD estimates (15 boundaries) | Derived by this code | Yes — `data/derived/boundary_estimates.csv` | CC-BY-4.0 |
| Pooled date × distance-bin panel (26,862 cells) | Derived by this code | Yes — `data/derived/panel_date_bin_pooled.csv` | CC-BY-4.0 |

### 2.3 Restricted data: what a replicator needs

The extract used in the paper contains, for each incident: date, time, incident category, and projected coordinates. It contains **no** direct personal identifiers — no names, no document numbers, no addresses. Coordinates are the location of the reported incident.

The extract can be requested from the Secretaría Distrital de Seguridad, Convivencia y Justicia through the Bogotá District Government's public information channel (*Bogotá te escucha*, Sistema Distrital para la Gestión de Peticiones Ciudadanas). `docs/data_access.md` reproduces the text of the request that produced the extract used here, the variable-level specification of the resulting file, and a checksum against which a newly obtained extract can be verified.

Typical turnaround at the time of writing was 15 working days. The author is not permitted to share the extract directly, including with referees or the journal's data editor; the derived files in `data/derived/` are provided for that purpose.

### 2.4 Human subjects and confidentiality

The records are de-identified administrative data on emergency incidents. No individual can be identified from the data or from any output in this package or the paper. No output is reported at a spatial resolution finer than a distance bin (50 m) aggregated over dates and boundaries; no incident location is published. See the Ethics statement in the manuscript.

---

## 3. Computational requirements

### 3.1 Software

- **R** 4.3.0 or later, with packages: `dplyr`, `tidyr`, `readr`, `purrr`, `tibble`, `ggplot2`, `sf` (≥ 1.0), `rdrobust` (≥ 2.2), `rddensity` (≥ 2.5)
- **Python** 3.10 or later, with `numpy`, `pandas`, `geopandas`, `matplotlib` (used only by `code/06_design_checks.py`)
- **GDAL/GEOS/PROJ**, required by `sf` and `geopandas`

Install the R dependencies with:

```r
install.packages(c("dplyr","tidyr","readr","purrr","tibble",
                   "ggplot2","sf","rdrobust","rddensity"))
```

`docs/session_info.txt` records the exact versions used to produce the reported results.

### 3.2 Hardware and runtime

Developed on macOS (Apple silicon), 16 GB RAM. The derived-data path runs in under two minutes. The full path takes 20–40 minutes, dominated by the point-in-polygon join and the per-date nearest-boundary distance computation in `02_build_panel.R`; peak memory is approximately 8 GB.

### 3.3 Randomness

No step depends on a random seed except the date block bootstrap in `code/99_extensions_rate_estimand.R`, which sets `set.seed(20240411)`. The randomization-inference test in `code/06_design_checks.py` is **not** Monte Carlo: it enumerates all 324 non-identity cyclic rotations of the study calendar exhaustively, so its *p*-value is exact and deterministic.

---

## 4. Instructions to replicators

1. Download the package and set the working directory to its root.
2. Open `code/_config.R`. On the derived-data path, change nothing. On the full path, set `RESTRICTED_INCIDENTS` to the path of your incident extract.
3. Run:

```bash
Rscript code/00_master.R
python3 code/06_design_checks.py
```

Every script prints its recomputed values next to the values reported in the paper, with a line marked `***` wherever the two diverge. A clean run prints no `***`.

Outputs are written to `output/tables/` and `output/figures/`. Nothing outside `output/` and `data/derived/` is modified.

---

## 5. Description of code

| File | Requires restricted data | Produces |
|---|---|---|
| `code/_config.R` | — | Paths, analysis constants, package check |
| `code/00_master.R` | — | Runs the pipeline |
| `code/01_build_adjacency.R` | No | Appendix A: adjacency pairs, degree statistics, shared boundary lengths |
| `code/02_build_panel.R` | **Yes** | Sample funnel (Section 3.4); boundary × date × bin panel |
| `code/03_estimate_boundary_rd.R` | **Yes** | Table 1 Panel A; bandwidths and effective samples |
| `code/04_aggregate_inference.R` | No | Table 1 Panels B–C; Section 5.2 equivalence and MDE; leave-one-out; Figures 3–4 |
| `code/05_robustness.R` | Partly | Table 2; Table 3; placebo cutoff estimates |
| `code/06_design_checks.py` | **Yes** | Table 4 rows 2–4; Figure 8 |
| `code/07_exhibits.R` | No | LaTeX table bodies for Tables 1 and C.1 |
| `code/08_python_reconstruction.py` | **Yes** | Independent Python rebuild of the whole chain; date-block bootstrap, exact Table 3 col (2), and the 08:00–08:00 re-timing check. Runs without R. |
| `code/99_extensions_rate_estimand.R` | **Yes** | Extensions not used for the reported estimates — see below |

### Note on `99_extensions_rate_estimand.R`

This script implements specifications discussed in the paper but **not** used for any reported estimate: the zero-inclusive, exposure-normalised Poisson rate estimand of Section 4.5.2; a date block bootstrap and Conley spatial-HAC standard errors for the pooled effect; a 2023 placebo-in-time; and a covariate-balance routine. Each requires an input the present version does not have — a residential population and *estrato* layer at 50-metre resolution, or georeferenced pre-rationing reports. The script is shipped so these specifications can be executed once those layers are obtained. **No number in the paper depends on it**, and it is deliberately excluded from `00_master.R`.

---

## 6. Description of data files

### `data/public/`

- **`hydraulic_polygons.gpkg`** — GeoPackage, EPSG:3116, 13 features. Twelve correspond to the EAAB service polygons in the rationing system; one additional geometry (identifier 4) is not assigned to any shift and contains no incidents. `01_build_adjacency.R` filters to the twelve using the crosswalk.
- **`rationing_calendar.csv`** — 325 rows. `date` (YYYY-MM-DD), `turno_group_id` (rationed shift, 1–9), `treated` (always 1; the file lists rationing days only). Days on which the programme was suspended — the June 2024 adjustment and 23 Dec 2024 – 6 Jan 2025 — are absent by construction.
- **`polygon_shift_crosswalk.csv`** — 12 rows: `polygon_id`, `shift_id`. Shift 5 = {6, 7}; Shift 8 = {1, 2, 3}; all others singletons. There is no polygon 4.

### `data/derived/`

- **`boundary_estimates.csv`** — 15 rows, one per estimation boundary. `boundary_id` (e.g. `3_11`), `poly1`, `poly2`, `shift1`, `shift2`, `n_crimes` (incidents within 2,000 m; sums to 20,729), `tau`, `se`, `pval`, `ci_lower`, `ci_upper`, `bandwidth` (MSE-optimal, metres), `n_eff_left`, `n_eff_right` (effective observations; sum to 4,451). **This is the file that reproduces the headline results.**
- **`boundary_estimates_with_counts.csv`** — the same estimates with `n_treated` and `n_control`, the within-bandwidth counts (summing to 5,951 and 8,309 respectively, i.e. 14,260 within-bandwidth observations). Used for Appendix Table C.1.
- **`panel_date_bin_pooled.csv`** — 26,862 rows: `date`, `bin_center` (metres, signed), `n_crimes`, `treated`, `side`. The pooled date × bin panel behind row 2 of Table 2.
- **`placebo_cutoffs_pooled.csv`** — 5 rows: pooled estimates at cutoff displacements of −500, −250, 0, +250 and +500 metres.
- **`unit_of_analysis_comparison.csv`** — 3 rows behind Table 2.
- **`adjacency_pairs.csv`, `adjacency_degree.csv`** — written by `01_build_adjacency.R`.

### A caution about superseded intermediate files

The author's working directory accumulated several intermediate result files over the project's history that are **not** the provenance of any reported number, including an earlier tabulated adjacency summary whose neighbour lists do not match the authoritative polygon layer. None of them is included here, and `01_build_adjacency.R` recomputes adjacency from the geometry at run time rather than reading any pre-tabulated list. Replicators working from the author's raw project directory rather than from this package should be aware of the distinction.

---

## 7. List of exhibits and their provenance

| Exhibit | Script | Data | Verified against paper |
|---|---|---|---|
| Table 1, Panel A | `03` → `07` | boundary estimates | Yes |
| Table 1, Panels B–C | `04` | boundary estimates | Yes — IVW −0.0438, SE 0.0567, Q 9.867, I² 0 |
| Table 2 | `05` | unit comparison | Yes |
| Table 3, cols (1)–(2) | `05` | boundary estimates | Yes — col (2) −0.0409, SE 0.0578, N 19,566, K 12 |
| Table 3, cols (3)–(4) | `05` (full path) | boundary × bin panel | Requires full path |
| Table 3, leave-one-out | `04` | boundary estimates | Yes — range [−0.062, −0.005] |
| Table 4, row 1 (density) | `05` (full path) | incident panel | Requires full path |
| Table 4, rows 2–4 | `06` | incident-level file | Requires full path |
| Table A.1 | — | crosswalk (public) | Yes |
| Table C.1 | `07` | boundary estimates | Yes |
| Figures 1, 2, A.1, A.2, C.1 | author's mapping code | polygons, calls | Maps and descriptive series |
| Figures 3, 4 | `04` | boundary estimates | Yes |
| Figure 5 | `04` | pooled estimate | Yes |
| Figures 6, 7 | `05` (full path) | boundary × bin panel | Requires full path |
| Figure 8 | `06` | incident-level file | Requires full path |

---

## 8. References

Calonico, S., M. D. Cattaneo and R. Titiunik (2014). "Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs." *Econometrica* 82(6), 2295–2326.

Cattaneo, M. D., M. Jansson and X. Ma (2018). "Manipulation Testing Based on Density Discontinuity." *Stata Journal* 18(1), 234–261.

Empresa de Acueducto y Alcantarillado de Bogotá (2024). *Resolución 291 de 2024*. Registro Distrital No. 7981, 10 April 2024.

Empresa de Acueducto y Alcantarillado de Bogotá (2024). *Resolución 980 de 2024*. 17 December 2024.

Secretaría Distrital de Seguridad, Convivencia y Justicia (2025). *Registros de incidentes — Línea de Emergencias 123, delitos contra la mujer, enero 2023 – mayo 2025* [restricted dataset]. Bogotá D.C.

---

## 9. Licence

Code in `code/` is released under the MIT Licence. Derived data in `data/derived/` is released under CC-BY-4.0. Public data in `data/public/` remains under the licence of its original publisher.
