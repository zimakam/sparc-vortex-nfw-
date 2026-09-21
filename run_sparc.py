#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SPARC VORTEX vs NFW: fitting script."""
import os, sys, json, math, time
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
RESULTS_DIR = os.path.join(HERE, "results")
os.makedirs(RESULTS_DIR, exist_ok=True)

MASS_FILE = os.path.join(DATA_DIR, "MassModels_Lelli2016c.mrt")
PROP_FILE = os.path.join(DATA_DIR, "SPARC_Lelli2016c.mrt")
OUT_CSV = os.path.join(RESULTS_DIR, "sparc_results.csv")
OUT_JSON = os.path.join(RESULTS_DIR, "summary.json")

UPS_B_RATIO = 1.4
A_FRACS = (0.0, 0.1, 0.2, 0.4, 0.8, 1.2)
RC_GRID = (0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 8.0)
UPS_GRID = (0.0, 0.4, 0.8, 1.0, 1.4, 2.0)

def parse_mass_models(path):
    g = defaultdict(list)
    with open(path) as f:
        dash = 0
        for line in f:
            if line.startswith("------"):
                dash += 1
                if dash >= 3: break
        for line in f:
            p = line.split()
            if len(p) < 8: continue
            try:
                g[p[0]].append(tuple(float(p[i]) for i in range(1, 8)))
            except ValueError:
                continue
    return g

def parse_properties(path):
    props = {}
    with open(path) as f:
        dash = 0
        for line in f:
            if line.startswith("------"):
                dash += 1
                if dash >= 3: break
        for line in f:
            s = line.strip()
            if not s or s.startswith("Note"): continue
            p = s.split()
            if len(p) < 18: continue
            if not all(c.isalnum() or c in "-_" for c in p[0]): continue
            try:
                props[p[0]] = {"Q": int(p[17])}
            except (ValueError, IndexError):
                continue
    return props

def v_vortex(r, A, rc):
    x = np.maximum(r / max(rc, 1e-9), 1e-9)
    f = np.maximum((x - np.arctan(x)) / x, 0.0)
    return A * np.sqrt(f)

def v_nfw(r, A_h, rs):
    x = np.maximum(r / max(rs, 1e-9), 1e-9)
    f = np.maximum(np.log(1.0 + x) - x / (1.0 + x), 0.0) / x
    return A_h * np.sqrt(f)

def model_v(mode, r, vg, vd, vb, params):
    A, rc, ups = params
    ub = ups * UPS_B_RATIO
    vbar2 = np.maximum(vg**2 + ups*vd**2 + ub*vb**2, 0.0)
    vh2 = (v_vortex(r, A, rc) if mode == "vortex" else v_nfw(r, A, rc))**2
    return np.sqrt(vbar2 + vh2)

def chi2(mode, r, vo, ve, vg, vd, vb, params):
    try:
        vm = model_v(mode, r, vg, vd, vb, params)
    except Exception:
        return 1e15
    if not np.all(np.isfinite(vm)):
        return 1e15
    err = np.maximum(ve, 0.02 * vo)
    c = float(np.sum(((vo - vm) / err)**2))
    return c if math.isfinite(c) else 1e15

def fit(mode, r, vo, ve, vg, vd, vb):
    vmax = float(np.max(vo))
    best_c = 1e15
    best = None
    for A in (a * vmax for a in A_FRACS):
        for rc in RC_GRID:
            for ups in UPS_GRID:
                c = chi2(mode, r, vo, ve, vg, vd, vb, (A, rc, ups))
                if c < best_c:
                    best_c = c
                    best = [float(A), float(rc), float(ups)]
    if best is None: return None
    if mode == "vortex":
        bounds = [(10.0, 500.0), (0.1, 20.0), (0.0, 2.0)]
    else:
        bounds = [(10.0, 500.0), (0.5, 100.0), (0.0, 2.0)]
    step = [max(abs(best[k]) * 0.3, 1e-3) for k in range(3)]
    for _ in range(200):
        improved = False
        for k in range(3):
            for s in (+1, -1):
                trial = list(best)
                trial[k] = min(max(trial[k] + s * step[k], bounds[k][0]), bounds[k][1])
                c = chi2(mode, r, vo, ve, vg, vd, vb, tuple(trial))
                if c < best_c - 1e-7:
                    best = trial
                    best_c = c
                    improved = True
        if not improved:
            step = [s * 0.5 for s in step]
            if max(step) < 1e-6:
                break
    N = len(r)
    dof = max(N - 3, 1)
    return {
        "mode": mode, "params": best,
        "chi2": best_c, "chi2_dof": best_c / dof,
        "bic": best_c + 3 * math.log(max(N, 2)),
        "aic": best_c + 6, "n_points": N,
    }

def main():
    print("=" * 68)
    print("SPARC VORTEX vs NFW")
    print("=" * 68)
    if not os.path.exists(MASS_FILE):
        print(f"ERROR: {MASS_FILE} not found")
        sys.exit(1)
    galaxies = parse_mass_models(MASS_FILE)
    props = parse_properties(PROP_FILE)
    print(f"Loaded {len(galaxies)} galaxies")
    results = []
    t0 = time.time()
    for name in sorted(galaxies):
        rows = galaxies[name]
        if len(rows) < 5: continue
        arr = np.array(rows)
        r, vo, ve, vg, vd, vb = (arr[:, i] for i in range(1, 7))
        rv = fit("vortex", r, vo, ve, vg, vd, vb)
        rn = fit("nfw", r, vo, ve, vg, vd, vb)
        if rv is None or rn is None: continue
        results.append({"name": name,
                        "Q": props.get(name, {}).get("Q", 0),
                        "n": len(rows), "vortex": rv, "nfw": rn})
    N = len(results)
    vc = np.array([r["vortex"]["chi2_dof"] for r in results])
    nc = np.array([r["nfw"]["chi2_dof"] for r in results])
    vb_sum = sum(r["vortex"]["bic"] for r in results)
    nb_sum = sum(r["nfw"]["bic"] for r in results)
    v_ok = int(np.sum(vc < 3.0))
    n_ok = int(np.sum(nc < 3.0))
    v_win = sum(1 for r in results if r["vortex"]["bic"] < r["nfw"]["bic"])
    print()
    print(f"  Galaxies:        {N}")
    print(f"  VORTEX:          {v_ok}/{N} ({100*v_ok/N:.1f}%)")
    print(f"    median chi2/dof: {np.median(vc):.3f}")
    print(f"    BIC sum:         {vb_sum:.0f}")
    print(f"  NFW:             {n_ok}/{N} ({100*n_ok/N:.1f}%)")
    print(f"    median chi2/dof: {np.median(nc):.3f}")
    print(f"    BIC sum:         {nb_sum:.0f}")
    print(f"  Delta BIC:       {nb_sum - vb_sum:+.0f} (in favour of VORTEX)")
    print(f"  Head-to-head:    VORTEX {v_win}, NFW {N - v_win}")
    with open(OUT_CSV, "w") as f:
        f.write("name,Q,n_points,"
                "v_chi2,v_chi2_dof,v_bic,v_A,v_rc,v_ups,"
                "n_chi2,n_chi2_dof,n_bic,n_A,n_rc,n_ups,delta_bic\n")
        for r in results:
            rv, rn = r["vortex"], r["nfw"]
            f.write(",".join([
                r["name"], str(r["Q"]), str(r["n"]),
                f"{rv['chi2']:.4f}", f"{rv['chi2_dof']:.4f}", f"{rv['bic']:.2f}",
                f"{rv['params'][0]:.4f}", f"{rv['params'][1]:.4f}", f"{rv['params'][2]:.4f}",
                f"{rn['chi2']:.4f}", f"{rn['chi2_dof']:.4f}", f"{rn['bic']:.2f}",
                f"{rn['params'][0]:.4f}", f"{rn['params'][1]:.4f}", f"{rn['params'][2]:.4f}",
                f"{nb_sum - vb_sum:.2f}"]) + "\n")
    summary = {
        "n_galaxies": N,
        "vortex": {"success_rate": v_ok/N, "n_ok": v_ok,
                   "median_chi2_dof": float(np.median(vc)),
                   "mean_chi2_dof": float(np.mean(vc)),
                   "bic_sum": float(vb_sum)},
        "nfw": {"success_rate": n_ok/N, "n_ok": n_ok,
                "median_chi2_dof": float(np.median(nc)),
                "mean_chi2_dof": float(np.mean(nc)),
                "bic_sum": float(nb_sum)},
        "delta_bic": float(nb_sum - vb_sum),
        "head_to_head": {"vortex_wins": v_win, "nfw_wins": N - v_win},
        "elapsed_seconds": round(time.time() - t0, 2)}
    with open(OUT_JSON, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\nSaved: {OUT_CSV}")
    print(f"Saved: {OUT_JSON}")
    print(f"Time: {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()
