#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Bootstrap 95% CI for Delta BIC."""
import os, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "results", "sparc_results.csv")) as f:
    rows = list(csv.DictReader(f))

vb = np.array([float(r["v_bic"]) for r in rows])
nb = np.array([float(r["n_bic"]) for r in rows])
delta = nb - vb
N = len(delta)

print(f"Galaxies: {N}")
print(f"Delta BIC (observed): {np.sum(delta):+.1f}")
print(f"Median per-galaxy:    {np.median(delta):+.2f}\n")

rng = np.random.default_rng(42)
B = 10000
idx = rng.integers(0, N, size=(B, N))
s = delta[idx].sum(axis=1)
m = np.median(delta[idx], axis=1)
w = (delta[idx] > 0).sum(axis=1)

lo, hi = np.percentile(s, [2.5, 97.5])
print(f"Delta BIC total:  {np.sum(delta):+.0f}")
print(f"  95% CI:         [{lo:+.0f}, {hi:+.0f}]\n")
lo, hi = np.percentile(m, [2.5, 97.5])
print(f"Delta BIC median: {np.median(delta):+.2f}")
print(f"  95% CI:         [{lo:+.2f}, {hi:+.2f}]\n")
lo, hi = np.percentile(w, [2.5, 97.5])
print(f"Wins VORTEX: {int((delta > 0).sum())}/{N}")
print(f"  95% CI:  [{lo:.0f}, {hi:.0f}]\n")
print(f"P(Delta BIC > 0): {np.mean(s > 0):.4f}")
