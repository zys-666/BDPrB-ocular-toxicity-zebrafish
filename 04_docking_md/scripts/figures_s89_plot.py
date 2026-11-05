# -*- coding: utf-8 -*-
"""Build Supplementary Figures S8 and S9."""
import numpy as np, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from PIL import Image

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
plt.rcParams["axes.linewidth"] = 1.2

OUT = r"D:\user data\Desktop\斑马鱼\14.补充材料"
BASE = r"D:\user data\Desktop\斑马鱼"
C = {"apo": "0.55", "BDPrB": "#82C0EA", "KyA33": "#F5C157"}


def style(ax):
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_linewidth(1.3)
    ax.tick_params(width=1.3, length=4, labelsize=9)


def load_xvg(p):
    return np.loadtxt(p, comments=("#", "@"))


# ======================= Figure S8 =======================
d = np.load(os.path.join(BASE, "5.dockMD", "figures", "s8_plddt.npz"))
res, pl = d["res"], d["plddt"]

fig = plt.figure(figsize=(7.4, 3.4))
gs = GridSpec(1, 2, width_ratios=[1.5, 1.0], wspace=0.22)
ax = fig.add_subplot(gs[0, 0])

ax.axhspan(0, 50, color="#F6C6C0", alpha=0.55, lw=0, label="very low (<50)")
ax.axhspan(50, 70, color="#FBE3B8", alpha=0.55, lw=0, label="low (50-70)")
ax.axhspan(70, 90, color="#D8ECF7", alpha=0.65, lw=0, label="confident (70-90)")
ax.axhspan(90, 100, color="#DDEEDC", alpha=0.65, lw=0, label="very high (>90)")
ax.axvspan(0.5, 25.5, color="white", alpha=0.0)
ax.plot(res, pl, lw=1.0, color="#2b2b2b")
ax.axvline(25.5, color="#C0392B", ls="--", lw=1.0)
ax.axvspan(0.5, 25.5, color="#E4574F", alpha=0.12, lw=0)
ax.set_xlim(0, 520); ax.set_ylim(0, 100)
ax.set_xlabel("Residue number", fontsize=11, fontweight="bold")
ax.set_ylabel("pLDDT", fontsize=11, fontweight="bold")
ax.text(280, 30, "mean pLDDT 98.2 (residues 26-518)", ha="center", fontsize=8, color="0.25")
ax.annotate("N-terminal 1-25 removed\n(mean pLDDT 41.0)", xy=(22, 42), xytext=(80, 22),
            fontsize=8, color="#8E2B22",
            arrowprops=dict(arrowstyle="->", color="#8E2B22", lw=0.9))
ax.legend(fontsize=6.5, frameon=False, loc="lower right", ncol=2, handlelength=1.2)
style(ax)

ax2 = fig.add_subplot(gs[0, 1])
img = Image.open(os.path.join(BASE, "5.dockMD", "figures", "s8_structure.png")).convert("RGB")
import numpy as _np
_a = _np.array(img.convert("L"))
_ys, _xs = _np.where(_a < 250)
_pad = 25
img = img.crop((max(0, _xs.min() - _pad), max(0, _ys.min() - _pad),
                min(img.width, _xs.max() + _pad), min(img.height, _ys.max() + _pad)))
ax2.imshow(img)
ax2.axis("off")
ax2.text(0.5, -0.02, "N-terminal 1-25 (red) reaches the catalytic cavity;\nCys320 shown as orange sphere",
         transform=ax2.transAxes, fontsize=7.5, color="0.2", ha="center", va="top")

fig.text(0.012, 0.95, "A", fontsize=13, fontweight="bold")
fig.text(0.635, 0.95, "B", fontsize=13, fontweight="bold")
plt.subplots_adjust(left=0.08, right=0.99, top=0.90, bottom=0.16)

for ext in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(OUT, "FigureS8.%s" % ext), dpi=600)
plt.close()
print("S8 saved")

# ======================= Figure S9 =======================
fig = plt.figure(figsize=(7.4, 6.6))
gs = GridSpec(3, 3, height_ratios=[1.0, 0.9, 0.9], hspace=0.55, wspace=0.42)

# --- A: H-bond occupancy ---
axA = fig.add_subplot(gs[0, 0])
pairs = ["L-N...\nASN475", "THR146...\nL-O", "CYS473...\nL-N2"]
vals = [3.0, 1.5, 0.6]
axA.bar(range(3), vals, color=C["KyA33"], edgecolor="black", lw=0.8, width=0.62)
axA.set_xticks(range(3)); axA.set_xticklabels(pairs, fontsize=6.5)
axA.set_ylabel("Occupancy (%)", fontsize=9, fontweight="bold")
axA.set_ylim(0, 5)
axA.set_title("Ligand-protein H-bonds", fontsize=9.5, fontweight="bold")
axA.text(0.5, 0.72, "BDPrB: 0%\n(no donor; ether O only)", transform=axA.transAxes,
         ha="center", fontsize=7, color="#1f5f8b")
style(axA)

# --- B: Rg ---
axB = fig.add_subplot(gs[0, 1])
for k, p in [("apo", "controls/apo/gyrate.xvg"), ("BDPrB", "MD/04_analysis/gyrate.xvg"),
             ("KyA33", "controls/kya33/gyrate.xvg")]:
    d = load_xvg(os.path.join(BASE, "5.dockMD", p))
    axB.plot(d[:, 0] / 1000.0, d[:, 1] * 10, lw=0.7, color=C[k], label=k)
axB.set_xlabel("Time (ns)", fontsize=9, fontweight="bold")
axB.set_ylabel("Rg (Å)", fontsize=9, fontweight="bold")
axB.set_title("Radius of gyration", fontsize=9.5, fontweight="bold")
axB.legend(fontsize=6.5, frameon=False, ncol=1, loc="lower right")
style(axB)

# --- C: SASA ---
axC = fig.add_subplot(gs[0, 2])
for k, p in [("apo", "controls/apo/sasa_prot.xvg"), ("BDPrB", "MD/03_production/sasa_prot.xvg"),
             ("KyA33", "controls/kya33/sasa_prot.xvg")]:
    d = load_xvg(os.path.join(BASE, "5.dockMD", p))
    axC.plot(d[:, 0] / 1000.0, d[:, 1], lw=0.7, color=C[k], label=k)
axC.set_xlabel("Time (ns)", fontsize=9, fontweight="bold")
axC.set_ylabel("SASA (nm$^2$)", fontsize=9, fontweight="bold")
axC.set_title("Solvent-accessible surface", fontsize=9.5, fontweight="bold")
axC.legend(fontsize=6.5, frameon=False, loc="lower right")
style(axC)

# --- D/E: contact heatmaps ---
for row, (tag, fn) in enumerate([("BDPrB", "s9_heat_BDPrB.npz"), ("KyA33", "s9_heat_KyA33.npz")]):
    ax = fig.add_subplot(gs[row + 1, :])
    d = np.load(os.path.join(BASE, "5.dockMD", "figures", fn))
    mat, rr = d["mat"], d["res"]
    k = np.ones(7) / 7.0
    sm = np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 1, mat).T
    im = ax.imshow(sm, aspect="auto", origin="lower", cmap="Blues",
                   extent=[rr.min(), rr.max(), 0, 100], vmin=0, vmax=1)
    ax.set_ylabel("Time (ns)", fontsize=8.5, fontweight="bold")
    ax.set_xlabel("Residue number", fontsize=9, fontweight="bold")
    ax.set_title("%s: per-residue contact occupancy" % tag, fontsize=9.5, fontweight="bold")
    ax.tick_params(labelsize=8)
    cb = plt.colorbar(im, ax=ax, pad=0.01, fraction=0.025)
    cb.set_label("contact", fontsize=7.5, labelpad=1)
    cb.ax.tick_params(labelsize=7)

fig.text(0.012, 0.965, "A", fontsize=13, fontweight="bold")
fig.text(0.352, 0.965, "B", fontsize=13, fontweight="bold")
fig.text(0.668, 0.965, "C", fontsize=13, fontweight="bold")
fig.text(0.012, 0.628, "D", fontsize=13, fontweight="bold")
fig.text(0.012, 0.335, "E", fontsize=13, fontweight="bold")
plt.subplots_adjust(left=0.09, right=0.90, top=0.93, bottom=0.07)

for ext in ("png", "svg", "pdf"):
    fig.savefig(os.path.join(OUT, "FigureS9.%s" % ext), dpi=600)
plt.close()
print("S9 saved")
