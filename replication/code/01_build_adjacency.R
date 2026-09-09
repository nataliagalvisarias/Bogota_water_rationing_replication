## ---------------------------------------------------------------------------
## 01_build_adjacency.R
## Reproduces: Appendix A (adjacency structure), Appendix Table A.1 statistics,
##             and the degree statistics quoted in Section 3.3.
## Inputs:  data/public/hydraulic_polygons.gpkg
##          data/public/polygon_shift_crosswalk.csv
## Outputs: data/derived/adjacency_pairs.csv
##          data/derived/adjacency_degree.csv
##
## This step uses only public data and is always run.
##
## IMPORTANT: adjacency is recomputed from the polygon geometry at run time.
## Do not substitute a pre-tabulated adjacency list: earlier vintages of this
## project contained a stale summary file whose neighbour lists do not match
## the authoritative layer, and the boundary set of Table 1 is not recoverable
## from it.
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(sf))

poly <- suppressWarnings(st_read(file.path(DIR_PUBLIC, "hydraulic_polygons.gpkg"), quiet = TRUE))
## The layer as originally distributed declared XYZ geometry while its
## sub-geometries were XY, which GDAL reports as corrupt and which makes
## st_transform() fail outright. Dropping any Z/M dimension first makes the read
## robust to either vintage of the layer.
poly <- st_zm(poly, drop = TRUE, what = "ZM")
poly <- st_transform(poly, CRS_METRIC)
poly <- poly[!st_is_empty(poly), ]
if (any(!st_is_valid(poly))) poly <- st_make_valid(poly)

## The layer carries one geometry per EAAB service polygon. Identify the
## polygon-id column defensively rather than assuming a name.
id_col <- intersect(c("polygon_id", "OBJECTID", "objectid", "id"), names(poly))[1]
if (is.na(id_col)) stop("Could not find a polygon identifier column in the layer.")
poly$polygon_id <- as.integer(poly[[id_col]])
stopifnot(!any(duplicated(poly$polygon_id)))

xwalk <- read_csv(file.path(DIR_PUBLIC, "polygon_shift_crosswalk.csv"),
                  show_col_types = FALSE)

## Keep only the twelve polygons that are part of the rationing system. The
## source layer contains one additional geometry (identifier 4) that is not
## assigned to any shift and contains no incidents; Section 2.2, footnote 2.
poly <- poly[poly$polygon_id %in% xwalk$polygon_id, ]
stopifnot(nrow(poly) == 12L)

## --- Contiguity: shared boundary of positive length --------------------------
## st_relate with the DE-9IM pattern "F***1****" selects pairs whose interiors
## are disjoint but whose boundaries share a one-dimensional (line) intersection.
## This is exactly the criterion in Appendix A: shared edges count, single
## touching points do not.
touching <- st_relate(poly, poly, pattern = "F***1****", sparse = TRUE)

pairs <- do.call(rbind, lapply(seq_along(touching), function(i) {
  js <- touching[[i]]
  js <- js[js > i]                       # each unordered pair once
  if (!length(js)) return(NULL)
  data.frame(poly_a = poly$polygon_id[i], poly_b = poly$polygon_id[js])
}))

## Shared boundary length, in metres, for each pair.
pairs$shared_length_m <- vapply(seq_len(nrow(pairs)), function(k) {
  a <- poly[poly$polygon_id == pairs$poly_a[k], ]
  b <- poly[poly$polygon_id == pairs$poly_b[k], ]
  as.numeric(st_length(st_intersection(st_boundary(a), st_boundary(b))))[1]
}, numeric(1))

## Drop degenerate intersections below the GIS tolerance (Appendix A, footnote).
pairs <- pairs[!is.na(pairs$shared_length_m) &
                 pairs$shared_length_m > ADJACENCY_TOL_M, ]

pairs <- pairs |>
  left_join(rename(xwalk, poly_a = polygon_id, shift_a = shift_id), by = "poly_a") |>
  left_join(rename(xwalk, poly_b = polygon_id, shift_b = shift_id), by = "poly_b") |>
  mutate(cross_shift = shift_a != shift_b) |>
  arrange(poly_a, poly_b)

write_csv(pairs, file.path(DIR_DERIVED, "adjacency_pairs.csv"))

## --- Degree statistics -------------------------------------------------------
deg <- bind_rows(
  transmute(pairs, polygon_id = poly_a),
  transmute(pairs, polygon_id = poly_b)
) |> count(polygon_id, name = "degree") |> arrange(polygon_id)
deg <- xwalk |> left_join(deg, by = "polygon_id") |>
  mutate(degree = tidyr::replace_na(degree, 0L))
write_csv(deg, file.path(DIR_DERIVED, "adjacency_degree.csv"))

## --- Verification against the paper ------------------------------------------
cat("\n--- Appendix A / Section 3.3 verification ---\n")
chk <- function(label, got, reported) {
  flag <- if (isTRUE(all.equal(got, reported, tolerance = 5e-3))) "OK " else "***"
  cat(sprintf("%s %-42s recomputed = %-10s paper = %s\n",
              flag, label, format(got, digits = 4), format(reported, digits = 4)))
}
## Expected values are those reported in Section 3.3 and Appendix A, which state the
## STRICT shared-edge rule. The project's QGIS working file records five further pairs
## under a snapping tolerance; none contributes identifying variation. See Appendix A.
chk("Adjacent pairs",                nrow(pairs),                    15)
chk("Pairs crossing a shift boundary", sum(pairs$cross_shift),       15)
chk("Mean polygon degree",           round(mean(deg$degree), 1),     2.5)
chk("SD of polygon degree",          round(sd(deg$degree), 1),       1.9)
chk("Minimum degree",                min(deg$degree),                0)
chk("Maximum degree",                max(deg$degree),                6)
chk("Shortest shared boundary (m)",  round(min(pairs$shared_length_m)),  1294)
chk("Longest shared boundary (m)",   round(max(pairs$shared_length_m)), 16853)
cat("Polygon with maximum degree: ", deg$polygon_id[which.max(deg$degree)],
    " (paper: 8)\n", sep = "")
cat("Rows marked *** indicate a divergence from the published value.\n")
