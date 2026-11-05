# -*- coding: utf-8 -*-
"""MD figures (SVG) for apo / BDPrB / KyA33."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "stix"

C = {"apo": "0.55", "BDPrB": "#1f77b4", "KyA33": "#e07b39"}
POCKET = [138, 141, 142, 145, 146, 188, 191, 192, 195, 262,
          314, 319, 320, 321, 475, 477, 478, 479, 483]

f = {
    "apo":   dict(rmsd="../controls/apo/rmsd_bb.xvg", rg="../controls/apo/gyrate.xvg",
                  rmsf="../controls/apo/rmsf_all.xvg", sasa="../controls/apo/pocket_sasa.xvg"),
    "BDPrB": dict(rmsd="../MD/04_analysis/rmsd_bb.xvg", rg="../MD/04_analysis/gyrate.xvg",
                  lig="../MD/04_analysis/rmsd_lig.xvg", md="../MD/03_production/mindist.xvg",
                  rmsf="../MD/03_production/rmsf_all.xvg", sasa="../MD/03_production/pocket_sasa.xvg"),
    "KyA33": dict(rmsd="../controls/kya33/rmsd_bb.xvg", rg="../controls/kya33/gyrate.xvg",
                  lig="../controls/kya33/rmsd_lig.xvg", md="../controls/kya33/mindist.xvg",
                  rmsf="../controls/kya33/rmsf_all.xvg", sasa="../controls/kya33/pocket_sasa.xvg"),
}

def load(p):
    return np.loadtxt(p, comments=("#", "@"))

def style(ax):
    ax.grid(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)

# ---------------- Fig 1: RMSD + Rg ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
for k in ("apo", "BDPrB", "KyA33"):
    d = load(f[k]["rmsd"])
    ax[0].plot(d[:, 0], d[:, 1] * 10, lw=0.7, color=C[k], label=k)
    g = load(f[k]["rg"])
    ax[1].plot(g[:, 0] / 1000.0, g[:, 1] * 10, lw=0.7, color=C[k], label=k)
ax[0].set_xlabel("Time (ns)"); ax[0].set_ylabel("Backbone RMSD (A)")
ax[0].set_title("Protein backbone RMSD", fontsize=9)
ax[0].legend(fontsize=7, frameon=False, ncol=3)
ax[1].set_xlabel("Time (ns)"); ax[1].set_ylabel("Radius of gyration (A)")
ax[1].set_title("Radius of gyration", fontsize=9)
for a in ax: style(a)
plt.tight_layout()
plt.savefig("fig1_rmsd_rg.svg")
plt.savefig("fig1_rmsd_rg.png", dpi=300)
plt.close()

# ---------------- Fig 2: ligand RMSD + min distance ----------------
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.8))
for k in ("BDPrB", "KyA33"):
    d = load(f[k]["lig"])
    ax[0].plot(d[:, 0], d[:, 1] * 10, lw=0.7, color=C[k], label=k)
    m = load(f[k]["md"])
    ax[1].plot(m[:, 0] / 1000.0, m[:, 1] * 10, lw=0.7, color=C[k], label=k)
ax[0].set_xlabel("Time (ns)"); ax[0].set_ylabel("Ligand RMSD (A)")
ax[0].set_title("Ligand RMSD (fit on protein)", fontsize=9)
ax[0].legend(fontsize=7, frameon=False)
ax[1].set_xlabel("Time (ns)"); ax[1].set_ylabel("Min ligand-protein distance (A)")
ax[1].set_title("Ligand-pocket contact", fontsize=9)
ax[1].set_ylim(0, 6)
for a in ax: style(a)
plt.tight_layout()
plt.savefig("fig2_ligand.svg")
plt.savefig("fig2_ligand.png", dpi=300)
plt.close()

# ---------------- Fig 3: RMSF ----------------
fig, ax = plt.subplots(figsize=(7.0, 3.0))
for k in ("apo", "BDPrB", "KyA33"):
    d = load(f[k]["rmsf"])
    ax.plot(d[:, 0], d[:, 1] * 10, lw=0.8, color=C[k], label=k)
for lo, hi in [(136, 148), (186, 196), (260, 264), (312, 323), (473, 485)]:
    ax.axvspan(lo, hi, color="#A5B59A", alpha=0.35, lw=0)
ax.set_ylim(0, 6)
ax.text(505, 5.6, "C-terminal tail up to 17 A (cut)", fontsize=6, ha="right", color="0.4")
ax.set_xlabel("Residue number"); ax.set_ylabel("RMSF (A)")
ax.set_title("Per-residue RMSF (pocket residues shaded)", fontsize=9)
ax.legend(fontsize=7, frameon=False, ncol=3)
ax.set_xlim(26, 518)
style(ax)
plt.tight_layout()
plt.savefig("fig3_rmsf.svg")
plt.savefig("fig3_rmsf.png", dpi=300)
plt.close()

# ---------------- Fig 4: pocket RMSF / SASA bars ----------------
pock_rmsf = {}
pock_sasa = {}
for k in ("apo", "BDPrB", "KyA33"):
    d = load(f[k]["rmsf"])
    res = d[:, 0].astype(int); v = d[:, 1] * 10
    pock_rmsf[k] = np.mean([v[res == x][0] for x in POCKET if (res == x).any()])
    s = load(f[k]["sasa"])
    pock_sasa[k] = s[:, 1].mean()

fig, ax = plt.subplots(1, 2, figsize=(6.4, 2.8))
ks = ["apo", "BDPrB", "KyA33"]
ax[0].bar(ks, [pock_rmsf[k] for k in ks], color=[C[k] for k in ks], width=0.6)
ax[0].set_ylabel("Pocket RMSF (A)"); ax[0].set_title("Pocket flexibility", fontsize=9)
ax[0].set_ylim(0, 1.6)
ax[1].bar(ks, [pock_sasa[k] for k in ks], color=[C[k] for k in ks], width=0.6)
ax[1].set_ylabel("Pocket SASA (nm²)"); ax[1].set_title("Pocket solvent exposure", fontsize=9)
ax[1].set_ylim(0, 32)
for a in ax: style(a)
plt.tight_layout()
plt.savefig("fig4_pocket_bars.svg")
plt.savefig("fig4_pocket_bars.png", dpi=300)
plt.close()

print("RMSF  ", {k: round(v, 2) for k, v in pock_rmsf.items()})
print("SASA  ", {k: round(v, 2) for k, v in pock_sasa.items()})
print("done")
