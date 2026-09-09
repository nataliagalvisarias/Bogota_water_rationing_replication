## ---------------------------------------------------------------------------
## 02_build_panel.R
## Reproduces: the sample-construction funnel of Section 3.4 and the
##             boundary x date x distance-bin estimation panel.
## Inputs:  RESTRICTED_INCIDENTS (see docs/data_access.md)
##          data/public/hydraulic_polygons.gpkg
##          data/public/rationing_calendar.csv
##          data/derived/adjacency_pairs.csv   (from step 01)
## Outputs: data/derived/panel_boundary_date_bin.csv
##          output/tables/sample_funnel.csv
##
## RUNS ONLY ON THE FULL PATH.  Without the restricted extract this script is
## skipped by 00_master.R and the shipped derived files are used instead.
## ---------------------------------------------------------------------------

if (is.na(RESTRICTED_INCIDENTS)) {
  stop("02_build_panel.R requires RESTRICTED_INCIDENTS to be set in _config.R.")
}
suppressPackageStartupMessages(library(sf))

## --- 1. Incidents ------------------------------------------------------------
## Expected columns in the extract (see docs/variable_definitions.md):
##   incident_id, date (YYYY-MM-DD), category, x, y
inc <- read_csv(RESTRICTED_INCIDENTS, show_col_types = FALSE)
required <- c("date", "category", "x", "y")
missing  <- setdiff(required, names(inc))
if (length(missing)) stop("Extract is missing columns: ", paste(missing, collapse = ", "))

n_raw <- nrow(inc)

inc <- inc |>
  mutate(date = as.Date(date)) |>
  filter(date >= STUDY_START, date <= STUDY_END,
         grepl("maltrato.*mujer", category, ignore.case = TRUE))
n_dv_window <- nrow(inc)

inc_sf <- st_as_sf(inc, coords = c("x", "y"), crs = CRS_METRIC, remove = FALSE)

## --- 2. Assign each incident to a hydraulic polygon --------------------------
poly  <- st_transform(st_read(file.path(DIR_PUBLIC, "hydraulic_polygons.gpkg"),
                              quiet = TRUE), CRS_METRIC)
xwalk <- read_csv(file.path(DIR_PUBLIC, "polygon_shift_crosswalk.csv"),
                  show_col_types = FALSE)
id_col <- intersect(c("polygon_id", "OBJECTID", "objectid", "id"), names(poly))[1]
poly$polygon_id <- as.integer(poly[[id_col]])
poly <- poly[poly$polygon_id %in% xwalk$polygon_id, ]

inc_sf <- st_join(inc_sf, poly["polygon_id"], join = st_within)
unmatched <- sum(is.na(inc_sf$polygon_id))
if (unmatched > 0) {
  warning(unmatched, " incidents fell outside every service polygon and were dropped. ",
          "If this exceeds a handful, check the coordinate reference system.")
  inc_sf <- inc_sf[!is.na(inc_sf$polygon_id), ]
}
inc_sf <- left_join(inc_sf, xwalk, by = "polygon_id")

## --- 3. Treatment status: equation (1) ---------------------------------------
cal <- read_csv(file.path(DIR_PUBLIC, "rationing_calendar.csv"), show_col_types = FALSE) |>
  transmute(date = as.Date(date), treated_shift = as.integer(turno_group_id))
inc_sf <- inc_sf |>
  inner_join(cal, by = "date") |>
  mutate(treated = as.integer(shift_id == treated_shift))

## --- 4. Adjacency screen: Section 3.4 step 1, equation (A.1) -----------------
pairs <- read_csv(file.path(DIR_DERIVED, "adjacency_pairs.csv"), show_col_types = FALSE)
## Undirected pair list -> directed neighbour lookup.
nbr <- bind_rows(
  transmute(pairs, polygon_id = poly_a, neighbour = poly_b),
  transmute(pairs, polygon_id = poly_b, neighbour = poly_a)
) |> left_join(rename(xwalk, neighbour = polygon_id, nbr_shift = shift_id),
               by = "neighbour")

## A control polygon is admissible on date t iff it is untreated and borders at
## least one polygon of the treated shift.
admissible <- nbr |>
  inner_join(cal, by = character(), relationship = "many-to-many") |>
  filter(nbr_shift == treated_shift) |>
  distinct(date, polygon_id)

inc_sf <- inc_sf |>
  left_join(mutate(admissible, admissible_control = TRUE),
            by = c("date", "polygon_id")) |>
  filter(treated == 1L | isTRUE(admissible_control) | !is.na(admissible_control))
n_after_adjacency <- nrow(inc_sf)

## --- 5. Nearest active boundary and signed running variable ------------------
## Active boundary set on date t: equation (7). Boundary geometry is the shared
## edge of the two polygons.
boundary_geom <- function(a, b) {
  st_intersection(st_boundary(poly[poly$polygon_id == a, ]),
                  st_boundary(poly[poly$polygon_id == b, ]))
}
bgeom <- lapply(seq_len(nrow(pairs)), function(k)
  boundary_geom(pairs$poly_a[k], pairs$poly_b[k]))
names(bgeom) <- paste(pairs$poly_a, pairs$poly_b, sep = "_")

assign_boundary <- function(day_df, day) {
  ts <- cal$treated_shift[cal$date == day]
  active <- which((pairs$shift_a == ts) != (pairs$shift_b == ts))   # exactly one side treated
  if (!length(active)) return(NULL)
  d <- vapply(active, function(k)
    as.numeric(st_distance(day_df, bgeom[[k]])), numeric(nrow(day_df)))
  d <- matrix(d, nrow = nrow(day_df))
  j <- max.col(-d, ties.method = "first")     # nearest active boundary, equation (4)
  day_df$boundary_id <- names(bgeom)[active[j]]
  day_df$dist_m      <- d[cbind(seq_len(nrow(d)), j)]
  day_df
}
panel_pts <- bind_rows(lapply(split(inc_sf, inc_sf$date), function(g)
  assign_boundary(g, unique(g$date))))

## --- 6. Distance screen and binning: Section 3.4 step 2, equation (10) -------
panel_pts <- panel_pts |> filter(dist_m <= MAX_DIST_M)
n_final <- nrow(panel_pts)

panel_pts <- panel_pts |>
  mutate(running   = ifelse(treated == 1L, -dist_m, dist_m),
         bin_index = floor(running / BIN_WIDTH_M + 0.5),
         bin_center = bin_index * BIN_WIDTH_M)

## Cell counts. Note: the reported specification counts incidents per occupied
## boundary x date x bin cell. The zero-inclusive variant that lays down the
## full bin grid is implemented in 99_extensions_rate_estimand.R; see
## Section 4.5.2 of the paper for why it is reported as an extension.
panel <- panel_pts |>
  st_drop_geometry() |>
  count(boundary_id, date, bin_center, treated, name = "n_crimes")

write_csv(panel, file.path(DIR_DERIVED, "panel_boundary_date_bin.csv"))

## --- 7. Sample funnel --------------------------------------------------------
funnel <- tibble::tibble(
  step = c("Raw extract rows",
           "Domestic violence, within study window",
           "After adjacency screen (Section 3.4 step 1)",
           "After 2,000 m distance screen (Section 3.4 step 2)"),
  n    = c(n_raw, n_dv_window, n_after_adjacency, n_final),
  reported = c(NA, 138130L, 50832L, 20729L)
)
write_csv(funnel, file.path(DIR_OUT_TAB, "sample_funnel.csv"))
print(funnel)
cat("\nThe 'reported' column gives the figures quoted in Section 3.4 and ",
    "Appendix A. Any divergence should be resolved before circulating results.\n", sep = "")
