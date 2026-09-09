#!/usr/bin/env python3
"""
09_rebuild_balanced_panel.py
=============================================================================
Rebuilds the boundary x date x bin panel from the incident-level RD estimation
sample and re-estimates the design under four nested specifications, isolating
the effect of each curation decision:

    S0  published specification        duplicates retained; bins from round(x/50)*50,
                                       so the cell at the cutoff spans both sides;
                                       occupied cells only
    S1  + deduplicated                 202 exact duplicate incident records dropped
    S2  + sign-preserving bins         bins defined as sign(x)*(floor(|x|/50)*50+25)
                                       so that no cell spans the cutoff
    S3  + balanced zero-filled panel   complete boundary x date x bin grid within
                                       +/-2,000 m, zero-incidence cells filled

It then runs a date-block bootstrap on S3 and rebuilds the unit-of-analysis
comparison under a common +/-2,000 m screen.

ESTIMATOR held constant across all specifications: local linear, triangular
kernel, MSE-optimal bandwidth, nearest-neighbour variance, clustered by date,
masspoints='off'.

IMPORTANT -- READ BEFORE CITING ANY NUMBER FROM THIS SCRIPT
    This uses the Python port of rdrobust. R's rdrobust is not installable on
    the R version available in this environment. Run against the published
    panel (S0) the port returns a pooled estimate of -0.0492 against the
    published -0.0438, with bandwidths agreeing to about 8 m on average and
    boundary estimates differing by 0.039 on average. It therefore reproduces
    the design closely but NOT exactly. The DIFFERENCES between S0, S1, S2 and
    S3 are computed under one engine and are the reportable content; the LEVELS
    should be regenerated with R's rdrobust before being written into the paper.

USAGE
    python3 09_rebuild_balanced_panel.py --sample rd_estimation_sample_attributes.csv \
                                         --hours  hours.csv --outdir .
REQUIREMENTS  pip install rdrobust pandas numpy scipy
=============================================================================
"""
import argparse, warnings
import numpy as np, pandas as pd
from scipy import stats
from rdrobust import rdrobust
warnings.filterwarnings("ignore")

BIN_W, MAX_D, N_BOOT, SEED = 50, 2000, 500, 20240411
RD = dict(c=0, p=1, kernel="triangular", bwselect="mserd", vce="nn", masspoints="off")

def bin_round(x):   return np.round(x / BIN_W).astype(int) * BIN_W
def bin_signed(x):  return np.sign(x) * (np.floor(np.abs(x) / BIN_W) * BIN_W + BIN_W / 2)

def cells_occupied(d, binf):
    d = d.copy(); d["bin"] = binf(d.running_var.values)
    return d.groupby(["boundary_id", "date", "bin"], as_index=False).size().rename(columns={"size": "n"})

def cells_filled(d, binf):
    """Complete balanced grid within +/-MAX_D for every observed boundary-date."""
    occ  = cells_occupied(d, binf)
    grid = np.concatenate([np.arange(-MAX_D + BIN_W/2, 0, BIN_W), np.arange(BIN_W/2, MAX_D, BIN_W)])
    pairs = occ[["boundary_id", "date"]].drop_duplicates()
    out = pairs.merge(pd.DataFrame({"bin": grid}), how="cross") \
               .merge(occ, on=["boundary_id", "date", "bin"], how="left")
    out["n"] = out.n.fillna(0).astype(int)
    return out

def fit_all(cells, label):
    rows = []
    for b, g in cells.groupby("boundary_id"):
        try:
            f = rdrobust(y=g.n.values, x=g["bin"].values.astype(float),
                         cluster=g.date.astype("int64").values, **RD)
            rows.append(dict(spec=label, boundary_id=b, tau=float(f.coef.iloc[0, 0]),
                             se=float(f.se.iloc[2, 0]), pval=float(f.pv.iloc[2, 0]),
                             h=float(f.bws.iloc[0, 0]), n_cells=len(g),
                             n_eff=int(f.N_h[0] + f.N_h[1]), n_incidents=int(g.n.sum())))
        except Exception as e:
            rows.append(dict(spec=label, boundary_id=b, tau=np.nan, se=np.nan,
                             n_cells=len(g), n_incidents=int(g.n.sum())))
            print(f"  [warn] {label}/{b}: {e}")
    return pd.DataFrame(rows)

def pooled(r):
    r = r.dropna(subset=["se"]); w = 1 / r.se ** 2
    t = (w * r.tau).sum() / w.sum(); s = 1 / np.sqrt(w.sum())
    Q = (w * (r.tau - t) ** 2).sum(); df = len(r) - 1
    return dict(tau=t, se=s, pval=2 * (1 - stats.norm.cdf(abs(t / s))),
                ci_lo=t - 1.96 * s, ci_hi=t + 1.96 * s, k=len(r), Q=Q,
                Q_p=1 - stats.chi2.cdf(Q, df), I2=max(0, (Q - df) / Q) * 100,
                n_incidents=int(r.n_incidents.sum()), n_cells=int(r.n_cells.sum()),
                n_eff=int(r.n_eff.sum()))

def bounds(tau, se, baseline):
    """Equivalence bounds expressed against the control-side mean of the SAME panel.
    Normalising by a baseline measured on a different unit (e.g. incidents per
    boundary-day, when tau is per date x bin cell) makes the percentages
    incomparable across specifications."""
    mde, hi, d = 2.80 * se, tau + 1.96 * se, 0.10 * baseline
    return dict(control_baseline_per_cell=baseline, mde=mde,
                mde_pct=100 * mde / baseline, upper_ci_pct=100 * hi / baseline,
                tost_p_10pct=max(1 - stats.norm.cdf((tau + d) / se), stats.norm.cdf((tau - d) / se)),
                tightest_equiv_pct=100 * (abs(tau) + stats.norm.ppf(0.95) * se) / baseline)

def date_block_bootstrap(panel, est, B=N_BOOT, seed=SEED):
    """Resample whole dates, holding bandwidths and IVW weights at full-sample values."""
    rng = np.random.default_rng(seed)
    bs = list(est.boundary_id); H = dict(zip(est.boundary_id, est.h)); W = 1 / est.se.values ** 2
    def icept(x, y, w):
        X = np.column_stack([np.ones_like(x), x]); XtW = X.T * w
        try: return np.linalg.solve(XtW @ X, XtW @ y)[0]
        except np.linalg.LinAlgError: return np.nan
    def tau_h(x, y, h):
        w = np.clip(1 - np.abs(x) / h, 0, None); L, R = (x < 0) & (w > 0), (x > 0) & (w > 0)
        if L.sum() < 3 or R.sum() < 3: return np.nan
        return icept(x[R], y[R], w[R]) - icept(x[L], y[L], w[L])
    pre, idx = {}, {}
    for b in bs:
        g = panel[panel.boundary_id == b]
        pre[b] = (g["bin"].values.astype(float), g.n.values.astype(float))
        idx[b] = pd.Series(range(len(g))).groupby(g.date.values).apply(list).to_dict()
    dates = np.array(sorted(panel.date.unique())); draws = np.empty(B)
    for i in range(B):
        s = rng.choice(dates, size=len(dates), replace=True); t = np.full(len(bs), np.nan)
        for k, b in enumerate(bs):
            ii = [j for dt in s for j in idx[b].get(dt, ())]
            if len(ii) >= 20:
                ii = np.array(ii); t[k] = tau_h(pre[b][0][ii], pre[b][1][ii], H[b])
        ok = ~np.isnan(t); draws[i] = (W[ok] * t[ok]).sum() / W[ok].sum()
    draws = draws[~np.isnan(draws)]
    return draws, float(1 / np.sqrt(W.sum()))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, help="rd_estimation_sample_attributes.csv")
    ap.add_argument("--hours",  required=False, help="crime_id,date,hour_of_day (for diagnostics)")
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()

    v = pd.read_csv(a.sample, low_memory=False); v["date"] = pd.to_datetime(v.date)
    v = v[v.in_estimation_sample == 1]
    dedup = v[v.is_duplicate_record == 0]
    print(f"incidents: {len(v):,} published / {len(dedup):,} deduplicated "
          f"({int(v.is_duplicate_record.sum())} exact duplicates dropped)")

    specs = {"S0_published":             cells_occupied(v,     bin_round),
             "S1_dedup":                 cells_occupied(dedup, bin_round),
             "S2_dedup_signedbins":      cells_occupied(dedup, bin_signed),
             "S3_dedup_signed_zerofill": cells_filled(dedup,   bin_signed)}

    per_b, summary = [], []
    for k, c in specs.items():
        r = fit_all(c, k); per_b.append(r)
        p = pooled(r); base = c.loc[c["bin"] > 0, "n"].mean()
        p.update(bounds(p["tau"], p["se"], base)); p["spec"] = k
        p["pct_zero_cells"] = 100 * (c.n == 0).mean()
        summary.append(p)
        print(f"{k:26s} tau={p['tau']:+.4f} SE={p['se']:.4f} p={p['pval']:.3f} "
              f"cells={p['n_cells']:,} baseline/cell={base:.3f} MDE={p['mde_pct']:.1f}% "
              f"equiv={p['tightest_equiv_pct']:.1f}%")

    pd.concat(per_b, ignore_index=True).to_csv(f"{a.outdir}/reestimation_by_boundary.csv", index=False)
    cols = ["spec","tau","se","pval","ci_lo","ci_hi","k","Q","Q_p","I2","n_incidents","n_cells",
            "n_eff","control_baseline_per_cell","mde","mde_pct","upper_ci_pct","tost_p_10pct",
            "tightest_equiv_pct","pct_zero_cells"]
    pd.DataFrame(summary)[cols].to_csv(f"{a.outdir}/reestimation_pooled.csv", index=False)
    specs["S3_dedup_signed_zerofill"].to_csv(f"{a.outdir}/panel_balanced_zerofilled.csv", index=False)

    est3 = pd.concat(per_b, ignore_index=True)
    est3 = est3[(est3.spec == "S3_dedup_signed_zerofill")].dropna(subset=["se"])
    panel3 = specs["S3_dedup_signed_zerofill"]
    draws, se_meta = date_block_bootstrap(panel3, est3)
    base3 = panel3.loc[panel3["bin"] > 0, "n"].mean()
    t3 = pooled(est3)["tau"]
    rows = []
    for lab, se in (("meta-analytic (independence assumed)", se_meta),
                    ("date-block bootstrap", float(draws.std(ddof=1)))):
        rows.append(dict(basis=lab, tau=t3, se=se, **bounds(t3, se, base3)))
    bt = pd.DataFrame(rows)
    bt.loc[bt.basis == "date-block bootstrap", "B"] = len(draws)
    bt.loc[bt.basis == "date-block bootstrap", "boot_ci_lower"] = np.percentile(draws, 2.5)
    bt.loc[bt.basis == "date-block bootstrap", "boot_ci_upper"] = np.percentile(draws, 97.5)
    bt.to_csv(f"{a.outdir}/bootstrap_zerofilled.csv", index=False)
    print(f"\nbootstrap (B={len(draws)}): SE {draws.std(ddof=1):.4f} vs meta-analytic {se_meta:.4f} "
          f"(ratio {draws.std(ddof=1)/se_meta:.3f})")
    print(f"\nWrote 4 CSVs to {a.outdir}/")

if __name__ == "__main__":
    main()
