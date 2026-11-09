# -*- coding: utf-8 -*-
"""Per-residue contribution bars for BDPrB and KyA33 (Graph2 style, system labelled)."""
import numpy as np, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
plt.rcParams["axes.linewidth"] = 1.2


def res_table(path):
    out = []
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) < 3 or p[0] == "Lig":
            continue
        nums = re.findall(r"[-+]?\d+\.\d+", line)
        if nums:
            out.append((p[1], float(nums[0])))
    return out


def draw(path, tag, fname):
    res = sorted(res_table(path), key=lambda t: t[1])[:15]
    names = ["P-" + r for r, _ in res]
    vals = [v for _, v in res]
    fig, ax = plt.subplots(figsize=(5.6, 3.7))
    col = "#82C0EA" if tag == "BDPrB" else "#F5C157"
    cols = [col] * len(vals)
    ax.bar(range(len(vals)), vals, color=cols, edgecolor="black", linewidth=0.8)
    ax.axhline(0, color="black", lw=1.0)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(names, rotation=90, fontsize=10)
    ax.set_ylabel("FEC (kJ/mol)", fontsize=13, fontweight="bold")
    ax.set_xlabel("Residues", fontsize=13, fontweight="bold")
    ax.set_title(tag, fontsize=13, fontweight="bold")
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.4)
    ax.tick_params(width=1.4, length=5, labelsize=11)
    plt.tight_layout()
    plt.savefig(fname + ".svg")
    plt.savefig(fname + ".png", dpi=300)
    plt.close()
    print(tag, [(r, round(v, 2)) for r, v in res[:5]])


draw("../MMPBSA/_pid~res_fixed.dat", "BDPrB", "mmpbsa_fig3_residues")
draw("../MMPBSA_kya33/_pid~res.dat", "KyA33", "mmpbsa_fig3_residues_KyA33")
