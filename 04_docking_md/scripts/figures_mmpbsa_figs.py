# -*- coding: utf-8 -*-
"""MM-PBSA figures in the style of Graph1/Graph2 (Times New Roman, no grid, black frames)."""
import numpy as np, re, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
plt.rcParams["axes.linewidth"] = 1.2

B = np.loadtxt("../MMPBSA/_pid~MMPBSA_fixed.dat", comments="#")     # time, dG, dG_DH, MM, COU, VDW, PB, SA
K = np.loadtxt("../MMPBSA_kya33/_pid~MMPBSA.dat", comments="#")
NAVY, RED = "#82C0EA", "#F5C157"

def frame(ax):
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.4)
    ax.tick_params(width=1.4, length=5, labelsize=11)

# ---------------- 1. dG distribution ----------------
fig, ax = plt.subplots(figsize=(4.6, 3.4))
bins = np.linspace(-190, -50, 30)
ax.hist(B[:, 2], bins=bins, color=NAVY, alpha=0.75, edgecolor="black", linewidth=0.8, label="BDPrB")
ax.hist(K[:, 2], bins=bins, color=RED, alpha=0.60, edgecolor="black", linewidth=0.8, label="KyA33")
ax.axvline(B[:, 2].mean(), color=NAVY, ls="--", lw=1.4)
ax.axvline(K[:, 2].mean(), color=RED, ls="--", lw=1.4)
ax.set_xlabel("Binding free energy (kJ/mol)", fontsize=13, fontweight="bold")
ax.set_ylabel("Frequency", fontsize=13, fontweight="bold")
ax.legend(frameon=False, fontsize=11)
frame(ax)
plt.tight_layout()
plt.savefig("mmpbsa_fig1_dG.svg")
plt.savefig("mmpbsa_fig1_dG.png", dpi=300)
plt.close()

# ---------------- 2. energy decomposition ----------------
labels = ["VDW", "COU", "PB", "SA", "dG"]
bv = [B[:, 5].mean(), B[:, 4].mean(), B[:, 6].mean(), B[:, 7].mean(), B[:, 2].mean()]
kv = [K[:, 5].mean(), K[:, 4].mean(), K[:, 6].mean(), K[:, 7].mean(), K[:, 2].mean()]
x = np.arange(len(labels)); w = 0.36
fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax.bar(x - w/2, bv, w, color=NAVY, edgecolor="black", linewidth=0.9, label="BDPrB")
ax.bar(x + w/2, kv, w, color=RED, edgecolor="black", linewidth=0.9, label="KyA33")
ax.axhline(0, color="black", lw=1.0)
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=12, fontweight="bold")
ax.set_ylabel("Contribution energy (kJ/mol)", fontsize=13, fontweight="bold")
ax.legend(frameon=False, fontsize=11)
frame(ax)
plt.tight_layout()
plt.savefig("mmpbsa_fig2_decomp.svg")
plt.savefig("mmpbsa_fig2_decomp.png", dpi=300)
plt.close()

print("done")
