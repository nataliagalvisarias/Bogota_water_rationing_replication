#!/usr/bin/env python3
"""
08_python_reconstruction.py
---------------------------------------------------------------------------
Independent Python reconstruction of the estimation pipeline, plus the three
checks that the R routine leaves as extensions:

    (a) date-block bootstrap of the pooled inverse-variance weighted estimate
    (b) exact recomputation of Table 3 column (2)
    (c) 08:00-08:00 treatment re-timing check

WHY THIS EXISTS.  The R pipeline (steps 01-07) is the provenance of the
reported estimates.  This script re-derives the whole chain -- polygon
assignment, adjacency, per-date boundary assignment, binning, local linear
RD -- from the incident microdata using only Python, as an independent check
and to run the three extensions above without an R installation.

VALIDATION STATUS.  Run against `code/crime_with_treatment_final.csv` this
script recovers the same fifteen estimation boundaries as the paper and
boundary-level estimates correlated 0.82 with the published ones (three match
to three decimals), but it does NOT reproduce the published sample exactly:
it finds 17,931 incidents within 2,000 m against the paper's 20,729, because
that CSV is a different vintage of the extract.  Results from this script are
therefore INDICATIVE.  Point it at the extract that produced Table 1 and the
numbers become directly reportable.

REQUIREMENTS
    pip install rdrobust shapely numpy pandas scipy

USAGE
    python3 08_python_reconstruction.py --incidents /path/to/extract.csv \
                                        --polygons  data/public/hydraulic_polygons.gpkg \
                                        --calendar  data/public/rationing_calendar.csv \
                                        --outdir    output/tables
---------------------------------------------------------------------------
"""
import argparse, re, sqlite3, sys, warnings
from itertools import combinations

import numpy as np
import pandas as pd
import shapely
from scipy import stats
from shapely import wkb

warnings.filterwarnings("ignore")

# --- analysis constants (mirror code/_config.R) ------------------------------
BIN_WIDTH_M     = 50
MAX_DIST_M      = 2000
ADJACENCY_TOL_M = 1.0
MIN_BOUNDARY_N  = 500
BASELINE_RATE   = 2.03
EXCLUDE_POLYGON = 4          # not part of the rationing system; see Section 2.2 fn 2
N_BOOT          = 500
SEED            = 20240411


# --- geopackage reader (avoids a geopandas/GDAL dependency) ------------------
def read_gpkg_polygons(path, table="new_shifts_keep"):
    """Return {polygon_id: {'geom': shapely geometry, 'shift': int}}."""
    con = sqlite3.connect(path)
    cols = [d[1] for d in con.execute(f"PRAGMA table_info('{table}')")]
    gi = cols.index("geom")
    oi = cols.index("OBJECTID")
    si = cols.index("racionam_1")            # holds "Shift N"

    def to_shape(blob):
        # GeoPackage BinaryHeader: magic(2) version(1) flags(1) srs_id(4) envelope
        env = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}[(blob[3] >> 1) & 0x07]
        return wkb.loads(blob[8 + env:])

    out = {}
    for r in con.execute(f"SELECT * FROM '{table}'"):
        pid = int(r[oi])
        out[pid] = {"geom": to_shape(r[gi]),
                    "shift": int(re.search(r"(\d+)", r[si]).group(1))}
    return out


def build_adjacency(P, keep):
    """Pairs of polygons sharing a boundary segment longer than the tolerance."""
    bnd = {p: P[p]["geom"].boundary for p in keep}
    pairs = []
    for a, b in combinations(keep, 2):
        if not P[a]["geom"].intersects(P[b]["geom"]):
            continue
        inter = bnd[a].intersection(bnd[b])
        if inter.length > ADJACENCY_TOL_M:
            pairs.append({"a": a, "b": b, "geom": inter, "length": inter.length})
    return pairs


def load_incidents(path):
    """
    The R export writes the point geometry as `c(x, y)`, which the CSV writer
    splits across two unquoted fields.  The header therefore has one fewer name
    than each data row, and naive parsing silently shifts every column by one.
    Read positionally with explicit names instead.
    """
    names = ["crime_id", "fecha", "hora", "nombre_dia", "anio", "mes", "tipo_cod",
             "tipo", "desc", "localidad", "barrio", "call_date", "polygon_id_old",
             "treated_old", "dist_old", "gx", "gy"]
    d = pd.read_csv(path, header=0, names=names, low_memory=False)
    d["x"] = pd.to_numeric(d.gx.str.replace("c(", "", regex=False))
    d["y"] = pd.to_numeric(d.gy.str.replace(")", "", regex=False))
    d["date"] = pd.to_datetime(d.call_date)
    d["hour"] = pd.to_numeric(d.hora, errors="coerce")
    if d.x.isna().any():
        sys.exit("Coordinate parsing failed; check the geometry column format.")
    return d


def assign_polygons(d, P, keep):
    pts = shapely.points(d.x.values, d.y.values)
    ids = np.array(keep)
    tree = shapely.STRtree(np.array([P[p]["geom"] for p in keep], dtype=object))
    hit = tree.query(pts, predicate="within")
    out = np.full(len(d), -1, dtype=int)
    out[hit[0]] = ids[hit[1]]
    return out


def build_panel(d, P, pairs, keep, cal, retime=False):
    """
    Assign each incident to its nearest ACTIVE boundary on its rationing day,
    apply the adjacency and distance screens, and bin the signed running
    variable.  With retime=True the rationing day runs 08:00-08:00, so an
    incident before 08:00 belongs to the previous day's suspension.
    """
    shift = {p: P[p]["shift"] for p in keep}
    d = d.copy()
    d["rday"] = (d.date - pd.to_timedelta((d.hour < 8).astype(int), unit="D")
                 if retime else d.date)
    d = d.merge(cal.rename(columns={"date": "rday"}), on="rday", how="inner")

    bid = lambda a, b: f"{min(a, b)}_{max(a, b)}"
    rows = []
    for _, g in d.groupby("rday", sort=True):
        ts = int(g.tshift.iloc[0])
        treated = {p for p in keep if shift[p] == ts}
        if not treated:
            continue
        active = [q for q in pairs
                  if (shift[q["a"]] == ts) != (shift[q["b"]] == ts)]
        if not active:
            continue
        controls = {q["a"] if shift[q["b"]] == ts else q["b"] for q in active}
        ok = g[g["poly"].isin(treated | controls)]
        if ok.empty:
            continue
        pts = shapely.points(ok.x.values, ok.y.values)
        D = np.column_stack([shapely.distance(pts, q["geom"]) for q in active])
        j = D.argmin(axis=1)
        rows.append(pd.DataFrame({
            "date": ok.rday.values,
            "treated": ok["poly"].isin(treated).astype(int).values,
            "boundary": [bid(active[k]["a"], active[k]["b"]) for k in j],
            "dist": D[np.arange(len(ok)), j]}))

    A = pd.concat(rows, ignore_index=True)
    A = A[A.dist <= MAX_DIST_M].copy()
    A["run"] = np.where(A.treated == 1, -A.dist, A.dist)
    A["bin"] = np.round(A.run / BIN_WIDTH_M).astype(int) * BIN_WIDTH_M
    return (A.groupby(["boundary", "date", "bin"], as_index=False)
              .size().rename(columns={"size": "n"}))


def fit_boundaries(cells, min_incidents=200):
    from rdrobust import rdrobust
    out = []
    for b, g in cells.groupby("boundary"):
        if g.n.sum() < min_incidents:
            continue                       # too sparse to support a local fit
        try:
            f = rdrobust(y=g.n.values, x=g.bin.values, c=0, p=1,
                         kernel="triangular", bwselect="mserd", vce="nn",
                         cluster=g.date.astype("int64").values)
            out.append({"boundary": b,
                        "tau": float(f.coef.iloc[0, 0]),
                        "se": float(f.se.iloc[2, 0]),
                        "pval": float(f.pv.iloc[2, 0]),
                        "h": float(f.bws.iloc[0, 0]),
                        "n_eff": int(f.N_h[0] + f.N_h[1]),
                        "n": int(g.n.sum())})
        except Exception as e:                       # noqa: BLE001
            print(f"  [warn] {b}: {e}")
    return pd.DataFrame(out)


def ivw(tau, se):
    w = 1.0 / np.asarray(se) ** 2
    t = float((w * np.asarray(tau)).sum() / w.sum())
    s = float(1.0 / np.sqrt(w.sum()))
    return t, s, 2 * (1 - stats.norm.cdf(abs(t / s))), t - 1.96 * s, t + 1.96 * s


# --- (a) date-block bootstrap -------------------------------------------------
def date_block_bootstrap(cells, est, B=N_BOOT, seed=SEED):
    """
    Resample whole DATES with replacement, preserving arbitrary within-date
    dependence across boundaries, and recompute the pooled estimate.

    Bandwidths and inverse-variance weights are held at their full-sample
    values.  Re-selecting the bandwidth inside each replicate would add
    selection noise that is not part of the sampling variability we want to
    measure, and is not what the meta-analytic standard error is being
    compared against.  The fixed-bandwidth weighted local linear fit below
    reproduces rdrobust's conventional coefficient exactly.
    """
    rng = np.random.default_rng(seed)
    bs = list(est.boundary)
    H = dict(zip(est.boundary, est.h))
    w0 = (1.0 / est.se.values ** 2)

    def intercept(x, y, w):
        X = np.column_stack([np.ones_like(x), x])
        XtW = X.T * w
        try:
            return np.linalg.solve(XtW @ X, XtW @ y)[0]
        except np.linalg.LinAlgError:
            return np.nan

    def tau_at_cutoff(x, y, h):
        w = np.clip(1 - np.abs(x) / h, 0, None)
        L, R = (x < 0) & (w > 0), (x >= 0) & (w > 0)
        if L.sum() < 3 or R.sum() < 3:
            return np.nan
        return intercept(x[R], y[R], w[R]) - intercept(x[L], y[L], w[L])

    pre, idx_by_date = {}, {}
    for b in bs:
        g = cells[cells.boundary == b]
        pre[b] = (g.bin.values.astype(float), g.n.values.astype(float))
        idx_by_date[b] = (pd.Series(range(len(g)))
                          .groupby(g.date.values).apply(list).to_dict())

    dates = np.array(sorted(cells.date.unique()))
    draws = np.empty(B)
    for r in range(B):
        samp = rng.choice(dates, size=len(dates), replace=True)
        taus = np.full(len(bs), np.nan)
        for k, b in enumerate(bs):
            idx = [i for dt in samp for i in idx_by_date[b].get(dt, ())]
            if len(idx) >= 10:
                idx = np.array(idx)
                taus[k] = tau_at_cutoff(pre[b][0][idx], pre[b][1][idx], H[b])
        m = ~np.isnan(taus)
        draws[r] = (w0[m] * taus[m]).sum() / w0[m].sum()

    draws = draws[~np.isnan(draws)]
    se_meta = float(1.0 / np.sqrt(w0.sum()))
    return {"se_bootstrap": float(draws.std(ddof=1)),
            "se_meta_analytic": se_meta,
            "ratio": float(draws.std(ddof=1) / se_meta),
            "ci_lower": float(np.percentile(draws, 2.5)),
            "ci_upper": float(np.percentile(draws, 97.5)),
            "B": int(len(draws))}


def equivalence_bounds(tau, se, baseline=BASELINE_RATE):
    mde = 2.80 * se
    hi = tau + 1.96 * se
    margin = 0.10 * baseline
    tost = max(1 - stats.norm.cdf((tau + margin) / se),
               stats.norm.cdf((tau - margin) / se))
    tight = abs(tau) + stats.norm.ppf(0.95) * se
    return {"mde_80pct": mde,
            "mde_pct_of_baseline": 100 * mde / baseline,
            "upper_ci_pct_of_baseline": 100 * hi / baseline,
            "tost_p_at_10pct_margin": tost,
            "tightest_equivalence_margin_pct": 100 * tight / baseline}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--incidents", required=True)
    ap.add_argument("--polygons", required=True)
    ap.add_argument("--calendar", required=True)
    ap.add_argument("--outdir", default="output/tables")
    ap.add_argument("--boot", type=int, default=N_BOOT)
    a = ap.parse_args()

    P = read_gpkg_polygons(a.polygons)
    keep = [p for p in sorted(P) if p != EXCLUDE_POLYGON]
    pairs = build_adjacency(P, keep)
    print(f"polygons kept: {len(keep)}   adjacency pairs: {len(pairs)}")

    cal = pd.read_csv(a.calendar)
    cal["date"] = pd.to_datetime(cal.date)
    cal["tshift"] = cal.turno_group_id.astype(int)
    cal = cal[["date", "tshift"]]
    print(f"rationing days: {len(cal)}  "
          f"({cal.date.min().date()} to {cal.date.max().date()})")

    inc = load_incidents(a.incidents)
    dv = inc[inc.desc.eq("MALTRATO A MUJER")].copy()
    dv = dv[(dv.date >= cal.date.min()) & (dv.date <= cal.date.max())]
    dv["poly"] = assign_polygons(dv, P, keep)
    dv = dv[dv["poly"] > 0]
    print(f"domestic-violence reports on rationing days: {len(dv):,}")

    # --- baseline (calendar-date treatment) ---------------------------------
    cells = build_panel(dv, P, pairs, keep, cal, retime=False)
    est = fit_boundaries(cells)
    t, s, p, lo, hi = ivw(est.tau, est.se)
    print(f"\nBASELINE  boundaries={len(est)}  incidents={int(cells.n.sum()):,}")
    print(f"  pooled IVW = {t:.4f}  SE = {s:.4f}  p = {p:.3f}  CI = [{lo:.4f}, {hi:.4f}]")

    # --- (a) bootstrap -------------------------------------------------------
    bb = date_block_bootstrap(cells, est, B=a.boot)
    print(f"\nDATE-BLOCK BOOTSTRAP (B={bb['B']})")
    print(f"  meta-analytic SE = {bb['se_meta_analytic']:.4f}")
    print(f"  bootstrap SE     = {bb['se_bootstrap']:.4f}   ratio = {bb['ratio']:.3f}")
    print(f"  percentile CI    = [{bb['ci_lower']:.4f}, {bb['ci_upper']:.4f}]")
    rows = [dict(basis="meta-analytic (independence assumed)", tau=t,
                 se=bb["se_meta_analytic"], **equivalence_bounds(t, bb["se_meta_analytic"])),
            dict(basis="date-block bootstrap", tau=t, se=bb["se_bootstrap"],
                 ci_lower=bb["ci_lower"], ci_upper=bb["ci_upper"], B=bb["B"],
                 **equivalence_bounds(t, bb["se_bootstrap"]))]
    pd.DataFrame(rows).to_csv(f"{a.outdir}/bootstrapped_se_results.csv", index=False)

    # --- (b) Table 3 column (2) ---------------------------------------------
    big = est[est.n >= MIN_BOUNDARY_N]
    t2, s2, p2, lo2, hi2 = ivw(big.tau, big.se)
    dropped = sorted(set(est.boundary) - set(big.boundary))
    print(f"\nTABLE 3 COL (2)  N_b >= {MIN_BOUNDARY_N}: {len(big)} boundaries, "
          f"{int(big.n.sum()):,} incidents; dropped {', '.join(dropped)}")
    print(f"  pooled IVW = {t2:.4f}  SE = {s2:.4f}  p = {p2:.3f}")
    pd.DataFrame([
        dict(specification="(1) baseline", tau=t, se=s, p_value=p, ci_lower=lo,
             ci_upper=hi, boundaries=len(est), incidents=int(cells.n.sum()),
             effective_obs=int(est.n_eff.sum()), excluded="none"),
        dict(specification=f"(2) N_b >= {MIN_BOUNDARY_N}", tau=t2, se=s2, p_value=p2,
             ci_lower=lo2, ci_upper=hi2, boundaries=len(big),
             incidents=int(big.n.sum()), effective_obs=int(big.n_eff.sum()),
             excluded=", ".join(dropped)),
    ]).to_csv(f"{a.outdir}/table3_col2_reconciled.csv", index=False)

    # --- (c) 08:00-08:00 re-timing ------------------------------------------
    cells_rt = build_panel(dv, P, pairs, keep, cal, retime=True)
    est_rt = fit_boundaries(cells_rt)
    tr, sr, pr, lor, hir = ivw(est_rt.tau, est_rt.se)
    share = 100 * dv.hour.lt(8).mean()
    print(f"\nRE-TIMED 08:00-08:00  ({share:.1f}% of incidents reassigned)")
    print(f"  pooled IVW = {tr:.4f}  SE = {sr:.4f}  p = {pr:.3f}  CI = [{lor:.4f}, {hir:.4f}]")
    print(f"  difference from calendar-date assignment: {tr - t:+.4f} "
          f"({abs(tr - t) / s:.2f} baseline standard errors)")
    head = pd.DataFrame([
        dict(level="POOLED", boundary="-", spec="calendar date", tau=t, se=s,
             p_value=p, ci_lower=lo, ci_upper=hi, incidents=int(cells.n.sum())),
        dict(level="POOLED", boundary="-", spec="re-timed 08:00-08:00", tau=tr, se=sr,
             p_value=pr, ci_lower=lor, ci_upper=hir, incidents=int(cells_rt.n.sum())),
        dict(level="POOLED", boundary="-", spec="difference", tau=tr - t),
    ])
    det = pd.concat([est.assign(level="boundary", spec="calendar date"),
                     est_rt.assign(level="boundary", spec="re-timed 08:00-08:00")],
                    ignore_index=True).rename(columns={"pval": "p_value", "n": "incidents"})
    pd.concat([head, det], ignore_index=True).to_csv(
        f"{a.outdir}/retimed_treatment_results.csv", index=False)

    print(f"\nWrote three CSVs to {a.outdir}/")


if __name__ == "__main__":
    main()
