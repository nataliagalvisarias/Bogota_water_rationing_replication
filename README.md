# Replication package

**Does Scheduled Urban Water Rationing Increase Reported Domestic Violence? Evidence from Bogotá’s Rotating Rationing Policy**

Natalia Galvis Arias, Global Development Institute, The University of Manchester.


- Repository: [Bogota_water_rationing_replication](https://github.com/nataliagalvisarias/Bogota_water_rationing_replication).
- Previous archived version `v1`: [10.5281/zenodo.22681258](https://doi.org/10.5281/zenodo.22681258), published 9 September 2026. It is **not** the archive of this candidate.
- Existing project concept DOI: [10.5281/zenodo.22681257](https://doi.org/10.5281/zenodo.22681257), covering the archive’s version series.
- Final version: `v2.0.0`. Archive DOI: [10.5281/zenodo.23700240](https://doi.org/10.5281/zenodo.23700240). The DOI must resolve after the Zenodo draft is published.
- Citation metadata: [CITATION.cff](CITATION.cff). The citation identifies the new version DOI; the older DOI below remains a historical reference.

## 1. Project overview and scope of replication

The paper studies Bogotá’s April 2024–April 2025 rotating water-rationing programme. Twelve hydraulic service polygons map to nine shifts. The design compares fixed 50-metre distance bands on opposite sides of twenty shared hydraulic boundaries. The primary model uses scheduled 08:00–08:00 operational water days, cell fixed effects and boundary-by-date fixed effects. The estimand is an intention-to-treat effect of scheduled assignment on **reported** domestic violence, rather than realised household water deprivation or underlying violence.

The shipped primary estimate is +0.0004 reports per identifying cell-day, with a 63-active-date moving-block bootstrap 95% interval of [−0.0130, +0.0135], or [−5.8%, +6.1%] of the comparison-side mean. These are rounded output values, not a claim of a statistically significant increase.

The default pipeline re-estimates statistical analyses from the shipped aggregate panels, regenerates nine statistical tables and seven statistical figures, and checks 47 displayed strings against the manuscript source. It does **not** independently rebuild the outcome panels from public microdata, rerun the original bandwidth selection, regenerate four manually composed maps, or verify every number in its contextual location. Access to restricted incidents is necessary for upstream reconstruction. See the limitations below.

## 2. Folder structure

```text
.
├── README.md                   scope, prerequisites and run instructions
├── CITATION.cff                candidate citation metadata
├── CODEBOOK.md                 data and output definitions
├── LICENSE                     MIT grant for code; data rights notice
├── .gitattributes              LF text endings and binary file rules
├── .gitignore                  raw data, session files and build exclusions
├── code/                       14 R scripts
├── data/
│   ├── public/                 schedule, source geometry and crosswalk
│   └── derived/                aggregate panels, areas, bandwidths and geography
├── figures/                    11 supplied manuscript exhibits
├── tables/                     9 generated LaTeX tables
├── output/                     43 shipped CSV results and diagnostics
├── manuscript/                 manuscript.tex and references.bib
└── docs/
    ├── DATA_ACCESS.md          restrictions, exact input schema, clearance status
    ├── DATA_SOURCES.md         source inventory and provenance gaps
    ├── SOFTWARE_REQUIREMENTS.md software versions and environment limits
    ├── VERIFICATION.md         verification evidence and its limitations
    ├── session_info.txt        supplied historical Linux environment record
    ├── audit_session_info_macos.txt current QC dependency session
    ├── FILE_MANIFEST.csv       paths, sizes, MD5, roles and review status
    └── checksums.md5           corresponding checksums
```

The manifest excludes its own two files and Git internals. It is an integrity inventory, not legal clearance or proof of statistical correctness. Internal review reports and temporary run logs are maintained outside this package.

## 3. Data sources, access and licensing

| Inputs | Source and role | Access and release status |
|---|---|---|
| Incident-level emergency calls, not included | Bogotá Secretaría Distrital de Seguridad, Convivencia y Justicia; C4/Línea 123; `MALTRATO`, with `MALTRATO A MUJER` for estimation | Restricted research extract. Request from custodian; approval is not guaranteed. |
| `data/public/rationing_calendar.csv` | Curated EAAB schedule: 306 assignments, 11 April 2024–11 April 2025; nine shifts | Included. Exact source documents, dates and reuse terms still need recording. |
| `data/public/hydraulic_polygons.gpkg` | EAAB hydraulic service polygons, EPSG:3116 | Included source geometry. Verify publisher’s licence and dataset/version identifier. |
| `data/public/polygon_shift_crosswalk.csv` | Mapping of 12 polygon IDs to nine shifts | Included curated input; document source and construction. |
| Five compressed count panels and monthly citywide series in `data/derived/` | Author-created aggregates of restricted calls | Author confirms these aggregates may be published. Small counts remain; individual calls are excluded. |
| Exposure areas, adjacency, repaired geometry and spatial weights | Geometry-only transformations | Included; retain EAAB provenance and applicable publisher terms. |
| `data/derived/bandwidths_mserd.csv` | Fifteen pinned historical bandwidths, plus median fallback for five boundaries | Included; original selection code/input/version not supplied. |

[DATA_SOURCES.md](docs/DATA_SOURCES.md) lists sources and unresolved provenance details. [CODEBOOK.md](CODEBOOK.md) describes variables, clocks, sample construction and outputs. Report–boundary placements are not unique incidents: a report near a junction can enter more than one fixed boundary panel.

No direct incident identifiers, incident coordinates or record-level times are intended for distribution. Aggregate cells have boundary, polygon, distance, date and, for segment panels, segment identifiers; these locate small counts indirectly. Geometry files contain public polygon coordinates. The permission statement relies on the author’s confirmation, rather than an inference based only on removal of coordinates.

**Licensing:** author-owned code is MIT; author-created derived aggregates are designated CC BY 4.0 in [LICENSE](LICENSE), consistent with the intended package licence. The author confirms aggregate publication is permitted. These grants do not cover the restricted incident extract or override municipal publishers’ terms. Record exact municipal reuse terms and source versions for archival completeness.

## 4. Software and computing requirements

- R **4.1 or later** is the syntax minimum; the supplied historical verification record specifies R 4.3.3 on Linux. A new Mac audit uses R 4.3.1. These are different environments.
- Default pipeline: `dplyr`, `readr`, `tidyr`, `sf`, `ggplot2`, `fixest`, `tibble`.
- Restricted rebuild additionally needs `lwgeom` and authorised native incident inputs.
- `sf` requires compatible GDAL, GEOS, PROJ and units libraries. CRAN binaries are preferable on supported macOS/Windows; Linux source builds may require system development libraries.
- `rdrobust` is not called by the shipped pipeline. Historical bandwidths are treated as fixed inputs.
- Optional manuscript compilation requires TeX Live or equivalent, `pdflatex`, BibTeX and the packages listed in [SOFTWARE_REQUIREMENTS.md](docs/SOFTWARE_REQUIREMENTS.md).
- No internet connection or credentials are needed during the default analytical run after dependencies are installed and the package is fully downloaded.

Install dependencies in R:

```r
install.packages(c("dplyr", "readr", "tidyr", "sf", "ggplot2", "fixest", "tibble"))
# Only for the authorised restricted-data path:
install.packages("lwgeom")
```

For a terminal installation, quoting that works in macOS/Linux shells is:

```bash
Rscript --vanilla -e 'install.packages(c("dplyr","readr","tidyr","sf","ggplot2","fixest","tibble"), repos="https://cloud.r-project.org")'
```

Windows users can run the R expressions in R/RStudio. Installation is not version-locked: `session_info.txt` records versions but is not an `renv.lock`, container, or exact restore mechanism. Do not describe the environment as frozen. Record `sessionInfo()` and `sf::sf_extSoftVersion()` for every verification run.

**Resources and runtime:** the supplied Linux note reports approximately eight minutes on one core, without a CPU model or peak-memory measurement. This is historical, not a guarantee for another computer. Script 10 performs 9,999 wild-bootstrap draws per weight family, multiple 1,000-draw block bootstraps and 999 permutations and is the main runtime cost. RAM and full-path runtime have not been benchmarked. The final QC report records the independent Mac outcome and timing. Keep a pristine copy and enough disk space for regenerated files.

## 5. Run the default aggregate pipeline

The default workflow uses the included aggregates; it requires no individual call records.

1. Copy or extract the complete package to a new working folder, retaining a pristine original. Make all cloud-backed files locally available. Use this folder as your working directory; running from its parent will fail.
2. Open a terminal in the package root (the folder containing this README), or set RStudio’s working directory there. Use the exact case of filenames. Paths within the pipeline use `file.path()` and relative paths.
3. Verify the original files **before** running anything that overwrites outputs:

```bash
Rscript --vanilla code/09_verify_manifest.R
```

A successful check prints `[09] verified ... files against the manifest`. Do not use `--refresh` to conceal an unexplained mismatch. The final review corrected locale-dependent path sorting so verification order is consistent across platforms.

4. Clear any unintended optional-layer configuration (`ADMIN_LAYERS`, `ESTRATO_LAYER`) and run in a fresh R process:

```bash
Rscript --vanilla code/00_master.R
```

For macOS/Linux, capture a run log **outside the package** to avoid changing its manifest membership:

```bash
Rscript --vanilla code/00_master.R > ../replication_run.log 2>&1
```

Windows PowerShell can use the same `Rscript --vanilla` invocation and redirect its output as appropriate. Default R expressions and relative file paths are portable; the environment-variable example in section 6 is shell-specific.

5. Require a successful process exit and these final messages:

```text
[08] all 47 displayed values found byte-for-byte in manuscript.tex
Replication complete. ...
```

Inspect warnings and compare regenerated statistical values/tables to the pristine copy. The 47-string check searches globally in the source; it does not prove correct labels, signs in context, completeness of reported results, or visual equivalence. Cross-platform floating-point noise may require numerical comparisons with declared tolerances; PDF metadata may differ even when rendered figures match.

The master overwrites regenerated CSVs, tables and statistical figures. It does not refresh the manifest. Only the release maintainer, after investigating all differences and finishing edits, should run:

```bash
Rscript --vanilla code/09_verify_manifest.R --refresh
Rscript --vanilla code/09_verify_manifest.R
```

Hash checks must be repeated on the exact downloaded release archive before submission.

## 6. Optional restricted-data rebuild

Authorised researchers must obtain an equivalent extract independently, retain it **outside** the package and confirm the field formats in [DATA_ACCESS.md](docs/DATA_ACCESS.md). A public open-data series is not established as equivalent to the restricted research extract.

In a macOS/Linux shell, from the package root:

```bash
INPUT_ROOT="/absolute/path/to/authorised/data" \
INCIDENTS_ALL="/absolute/path/to/authorised/data/raw/NUSE_NGA.dsv" \
Rscript --vanilla code/00_master.R --full
```

In Windows PowerShell:

```powershell
$env:INPUT_ROOT = "C:/authorised/data"
$env:INCIDENTS_ALL = "C:/authorised/data/raw/NUSE_NGA.dsv"
Rscript --vanilla code/00_master.R --full
```

`INPUT_ROOT` must be an existing directory and `INCIDENTS_ALL` an existing `.dsv` or supported `.gpkg`. Scripts 02 and 03 reconstruct the fixed-geography panels and areas; scripts 05 and 06 reconstruct treatment clocks and segment panels; script 11 adds restricted diagnostics. Outputs go below the package root. No raw record should be added to Git, archives or run logs.

The `.dsv` adapter expects semicolon separation; incident IDs; `FECHA` parsed with `%d-%b-%y`; integer hours `0`–`23`; exact Spanish classification labels; WGS84 longitude/latitude. Month abbreviations depend on `LC_TIME`. Inspect date parse failures, duplicate IDs and source coverage before accepting a restricted rebuild. Deduplication retains one row per incident ID, not merely exact duplicate rows. The restricted and external-layer paths were not executed in the final audit because their inputs were unavailable.

## 7. Execution order

Always prefer the master: later scripts reuse functions/objects defined earlier.

| Order | Script | Function |
|---|---|---|
| Pre-run | `09_verify_manifest.R` | Integrity check; not called by master. |
| 1 | `00a_prepare_geography.R` | Allocate narrow overlaps by ascending polygon ID; validate 20-edge membership under reversed priority. |
| 2 | `01_validate_adjacency.R` | Shared-edge adjacency, weights, neighbour list. |
| Full only | `02_read_restricted_incidents.R` | Read and filter incidents in memory; monthly descriptive count export. |
| Full only | `03_build_fixed_geography.R` | Fixed cell panels, areas and assignment diagnostics. |
| 3 | `04_placebo_checks.R` | Spatial and ancillary shift-label placebos. |
| 4 | `05_treatment_clocks.R` | Operational/calendar models, local boundary fits, block-bootstrap inference. |
| 5 | `06_robustness.R` | Dynamics, common bandwidth, segments, PPML and coverage. |
| 6 | `10_additional_analyses.R` | Wider dynamics, wild inference, block sensitivity, heterogeneity and primary-specification placebos. |
| Full only | `11_restricted_additional_analyses.R` | Additional restricted diagnostics. |
| Configured only | `12_external_layers.R` | Optional administrative overlaps and estrato analyses. |
| 7 | `07_manuscript_exhibits.R` | Regenerate statistical tables and figures. |
| 8 | `08_verify_manuscript_numbers.R` | 47 displayed-string presence checks. |

`00_master.R` orchestrates this order and checks the calendar and prerequisite inputs/packages. No optional output from scripts 11/12 is needed for the reported statistical exhibits. Shipped panels and the monthly series nevertheless depend upstream on restricted data.

## 8. Outputs and manuscript correspondence

| Exhibit or result | Source / regeneration |
|---|---|
| Figures 1, A1, A2, C2 | Supplied QGIS maps; **not regenerated**. No QGIS project, styles or exact export instructions supplied. |
| Figure 2 | Script 07; `data/derived/monthly_maltrato_citywide.csv` (broader citywide MALTRATO series). |
| Figures 3–5 | Script 07; boundary estimates, distance profile and primary timing outputs. |
| Figures 6–7 | Script 07; operational dense placebos and Estimator B shift permutations. |
| Figure C1 | Script 07; wider event-time output. |
| Tables 1–4, C1–C5 | Script 07; generated numerical CSVs. |
| Table A1 | Typeset inside manuscript; check against `polygon_shift_crosswalk.csv`. |
| Primary/secondary estimates | `output/primary_timing.csv`, `boundary_estimates_operational.csv`, `boundary_estimates_calendar.csv`. |
| Robustness and inference | `dynamic_*.csv`, `common500_*.csv`, `along_boundary_segments.csv`, `ppml.csv`, `ext_*.csv`. |
| Geography and construction | `adjacency_candidates_audit.csv`, `topology_*.csv`, `bandwidth_provenance.csv`, coverage outputs. |
| Manuscript check | `output/manuscript_number_audit.csv`; not a comprehensive statistical audit. |

Detailed variables and output families are in [CODEBOOK.md](CODEBOOK.md).

## 9. Optional manuscript compilation

Run from `manuscript/` after generating exhibits:

```bash
pdflatex -halt-on-error manuscript.tex
bibtex manuscript
pdflatex -halt-on-error manuscript.tex
pdflatex -halt-on-error manuscript.tex
```

Inspect undefined references/citations and the resulting PDF. The revised packaged source compiled successfully with the desktop editor’s built-in compiler during final QC. No new PDF was exported or visually certified. The historical 56-page claim was not re-counted; regenerate and review the submission PDF after all final edits. Do not run compilation before verifying the untouched manifest, since it creates build files. Before release, remove build products and verify again.

## 10. Replication notes and remaining limitations

- This is an aggregate-input replication package, with an optional restricted reconstruction path; it is not a fully public raw-to-results workflow.
- The author confirms permission to publish the included aggregates. A source-terms/reference note would strengthen the provenance record; no signed agreement has been independently examined in this review. Individual records remain restricted.
- Missing outcome dates are dropped rather than automatically treated as true zero-report dates. The operational clock has 304 usable dates; the calendar clock has 299 out of 306 schedule assignments. The distinction between source gaps and actual zero-report dates needs documented source coverage evidence.
- Fifteen bandwidths are pinned; five use their median, about 480.86 m. Original bandwidth-selection code and environment are absent.
- Maps are final static exhibits; editable map projects/export recipes are absent.
- Public-input provenance lacks exact download dates, catalogue/version identifiers and documented reuse terms.
- A recorded software session is supplied, but no complete dependency lockfile or container. Windows has not been independently tested.
- Manifest integrity and 47 literal string matches are useful checks with limited scope. Neither establishes confidentiality, identification validity or journal compliance.
- Existing archive `v1` describes an older paper title and 15 boundaries. Its metadata and files must not be represented as this 20-boundary candidate.

## 11. Citation and contact

Cite: Galvis Arias, Natalia. 2026. Replication Data and Code for: Does Scheduled Urban Water Rationing Increase Reported Domestic Violence? Evidence from Bogotá’s Rotating Rationing Policy. Version 2.0.0. Zenodo. [10.5281/zenodo.23700240](https://doi.org/10.5281/zenodo.23700240). The concept DOI can identify the project series but cannot identify the exact submitted files. The previous version DOI identifies only `v1`.

**Contact:** Natalia Galvis Arias, [natalia.galvisarias@postgrad.manchester.ac.uk](mailto:natalia.galvisarias@postgrad.manchester.ac.uk), Global Development Institute, The University of Manchester. [ORCID 0009-0004-1346-5138](https://orcid.org/0009-0004-1346-5138).

Report non-sensitive code issues through the [repository issue tracker](https://github.com/nataliagalvisarias/Bogota_water_rationing_replication/issues). Send access and confidentiality questions privately to the author or data custodian; do not post incident information in public issues.
