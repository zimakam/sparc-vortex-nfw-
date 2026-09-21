#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Generate 4 publication figures."""
import os, csv, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results")
with open(os.path.join(R, "sparc_results.csv")) as f:
    rows = list(csv.DictReader(f))
with open(os.path.join(R, "summary.json")) as f:
    summary = json.load(f)

vc = np.array([float(r["v_chi2_dof"]) for r in rows])
nc = np.array([float(r["n_chi2_dof"]) for r in rows])
vb = np.array([float(r["v_bic"]) for r in rows])
nb = np.array([float(r["n_bic"]) for r in rows])
delta = nb - vb
Q = np.array([int(r["Q"]) for r in rows])

fig, ax = plt.subplots(figsize=(8, 5))
bins = np.linspace(0, 8, 40)
ax.hist(vc, bins=bins, alpha=0.65, color="#2E86AB",
        label=f"VORTEX (median={np.median(vc):.2f})")
ax.hist(nc, bins=bins, alpha=0.65, color="#A23B72",
        label=f"NFW (median={np.median(nc):.2f})")
ax.axvline(1.0, ls="--", color="gray", alpha=0.6)
ax.axvline(3.0, ls=":", color="red", alpha=0.6)
ax.set_xlabel(r"$\chi^2/\mathrm{dof}$", fontsize=12)
ax.set_ylabel("Number of galaxies", fontsize=12)
ax.set_title(r"Reduced $\chi^2$ distribution: VORTEX vs NFW", fontsize=13)
ax.legend(loc="upper right"); ax.set_xlim(0, 8)
plt.tight_layout()
plt.savefig(os.path.join(R, "fig1_chi2_distribution.png"), dpi=150)
plt.close()

fig, ax = plt.subplots(figsize=(8, 5))
ax.hist(delta, bins=40, color="#4C9F70", edgecolor="black", linewidth=0.5)
ax.axvline(0, color="black", linewidth=1)
ax.axvline(np.sum(delta), color="blue", ls=":", label=f"total = {np.sum(delta):+.0f}")
ax.set_xlabel(r"$\Delta\mathrm{BIC}$", fontsize=11)
ax.set_ylabel("Number of galaxies", fontsize=12)
ax.set_title("Per-galaxy BIC difference", fontsize=13)
ax.legend(); ax.set_xlim(-100, 300)
plt.tight_layout()
plt.savefig(os.path.join(R, "fig2_bic_histogram.png"), dpi=150)
plt.close()

fig, ax = plt.subplots(figsize=(7, 7))
sc = ax.scatter(vc, nc, c=Q, cmap="viridis", s=30, alpha=0.7,
                edgecolor="black", linewidth=0.3)
lim = [0, min(10, max(vc.max(), nc.max()) * 1.05)]
ax.plot(lim, lim, "k--", alpha=0.5, label="1:1")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel(r"VORTEX $\chi^2/\mathrm{dof}$", fontsize=12)
ax.set_ylabel(r"NFW $\chi^2/\mathrm{dof}$", fontsize=12)
ax.set_title("Per-galaxy comparison", fontsize=13)
ax.legend(loc="upper left")
plt.colorbar(sc, ax=ax, label="Quality flag Q")
plt.tight_layout()
plt.savefig(os.path.join(R, "fig3_vortex_vs_nfw.png"), dpi=150)
plt.close()

fig, ax = plt.subplots(figsize=(7, 5))
succ = [100 * summary["vortex"]["success_rate"],
        100 * summary["nfw"]["success_rate"]]
bars = ax.bar(["VORTEX", "NFW"], succ, color=["#2E86AB", "#A23B72"],
              width=0.5, edgecolor="black")
for b, v in zip(bars, succ):
    ax.text(b.get_x() + b.get_width()/2, v + 1.5,
            f"{v:.1f}%", ha="center", fontweight="bold")
ax.set_ylabel(r"Success rate ($\chi^2/\mathrm{dof}<3$), %")
ax.set_title(f"Success rate over {summary['n_galaxies']} SPARC galaxies")
ax.set_ylim(0, 100); ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(R, "fig4_success_rate.png"), dpi=150)
plt.close()

print("Generated 4 figures in", R)
