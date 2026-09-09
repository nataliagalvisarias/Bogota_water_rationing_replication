# Variable definitions

## Restricted extract (not shipped)

Expected by `code/02_build_panel.R`. Column names may differ in a newly obtained
extract; rename to match before running, or adjust the `required` vector in that script.

| Column | Type | Definition |
|---|---|---|
| `incident_id` | string | Record identifier assigned by the 123 system. Not used in estimation; retained for de-duplication. |
| `date` | date (YYYY-MM-DD) | Calendar date of the incident report. |
| `time` | time (HH:MM:SS) | Time of the report. Not used in the baseline, which assigns treatment by calendar date; required for the 08:00–08:00 re-timing discussed in Section 4.1, footnote 3. |
| `category` | string | Incident classification. The main outcome filters on *maltrato a mujer* (domestic violence, female victim). The comparison outcome in Table 4 filters on sexual violence. |
| `x`, `y` | numeric | Projected coordinates, EPSG:3116 (MAGNA-SIRGAS / Colombia Bogotá zone), metres. If the extract supplies geographic coordinates (EPSG:4326) instead, reproject before use — distances in the running variable must be metric. |

## Derived: `boundary_estimates.csv`

| Column | Definition |
|---|---|
| `boundary_id` | Polygon pair, lower id first, e.g. `3_11`. |
| `poly1`, `poly2` | EAAB polygon identifiers of the two sides. |
| `shift1`, `shift2` | Rationing shifts to which those polygons belong. Always different. |
| `n_crimes` | Incidents assigned to this boundary within 2,000 m over the study period. Sums to 20,729. |
| `tau` | Boundary-specific RD estimate, incidents per date × 50 m bin. Conventional local-linear coefficient. |
| `se` | Robust bias-corrected standard error, nearest-neighbour variance, clustered by date. |
| `pval`, `ci_lower`, `ci_upper` | Robust bias-corrected inference. |
| `bandwidth` | MSE-optimal bandwidth, metres. |
| `n_eff_left`, `n_eff_right` | Effective observations on each side of the cutoff. Sum to 4,451 across boundaries. |

## Derived: `panel_date_bin_pooled.csv`

| Column | Definition |
|---|---|
| `date` | Calendar date. |
| `bin_center` | Signed distance-bin midpoint in metres. Negative = treated side, positive = control side (equation 6). |
| `n_crimes` | Incident count in the cell. |
| `treated` | 1 if the cell is on the treated side. |
| `side` | Human-readable label for `treated`. |

## Sign convention

Throughout, the running variable is negative on the **treated** side and positive on the
**control** side. A positive $\hat\tau_b$ therefore means *more* reported violence under
rationing. See equation (6) and Section 4.3.2.
