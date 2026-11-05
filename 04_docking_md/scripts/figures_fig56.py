# -*- coding: utf-8 -*-
"""Contact-frequency figure + snapshot strip."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "stix"
from PIL import Image, ImageDraw, ImageFont

CB = "../MD/03_production/contacts_BDPrB.csv"
CK = "../controls/kya33/contacts_KyA33.csv"

def load(p):
    d = np.loadtxt(p, comments="#")
    return d[:, 0].astype(int), d[:, 1]

r1, v1 = load(CB)
r2, v2 = load(CK)
m1 = dict(zip(r1.tolist(), v1)); m2 = dict(zip(r2.tolist(), v2))

POCKET = [(136, 148), (186, 196), (260, 264), (312, 323), (473, 485)]
C = {"BDPrB": "#1f77b4", "KyA33": "#e07b39"}

fig, ax = plt.subplots(1, 2, figsize=(7.4, 2.9), gridspec_kw={"width_ratios": [2.2, 1]})
for lo, hi in POCKET:
    ax[0].axvspan(lo, hi, color="#A5B59A", alpha=0.30, lw=0)
ax[0].step(r1, v1, where="mid", lw=0.8, color=C["BDPrB"], label="BDPrB")
ax[0].step(r2, v2, where="mid", lw=0.8, color=C["KyA33"], label="KyA33")
ax[0].set_xlim(26, 518); ax[0].set_ylim(0, 105)
ax[0].set_xlabel("Residue number"); ax[0].set_ylabel("Contact occupancy (%)")
ax[0].set_title("Per-residue ligand contact occupancy", fontsize=9)
ax[0].legend(fontsize=7, frameon=False, ncol=2)
ax[0].grid(False)

top = sorted(m1.items(), key=lambda kv: -kv[1])[:10]
labels = ["%d" % r for r, _ in top]
x = np.arange(len(top))
w = 0.4
ax[1].bar(x - w/2, [m1[r] for r, _ in top], w, color=C["BDPrB"], label="BDPrB")
ax[1].bar(x + w/2, [m2.get(r, 0) for r, _ in top], w, color=C["KyA33"], label="KyA33")
ax[1].set_xticks(x); ax[1].set_xticklabels(labels, rotation=60, fontsize=6)
ax[1].set_ylim(0, 110); ax[1].set_ylabel("Occupancy (%)")
ax[1].set_title("Top contact residues (BDPrB)", fontsize=9)
ax[1].legend(fontsize=7, frameon=False)
ax[1].grid(False)
for a in ax:
    a.spines["top"].set_visible(False); a.spines["right"].set_visible(False)
plt.tight_layout()
plt.savefig("fig5_contacts.png", dpi=300)
plt.savefig("fig5_contacts.svg")
plt.close()

# ---- snapshot strip ----
files = [("snap/pocket_0.png", "0 ns"), ("snap/pocket_25.png", "25 ns"),
         ("snap/pocket_50.png", "50 ns"), ("snap/pocket_75.png", "75 ns"),
         ("snap/pocket_100.png", "100 ns")]
imgs = [Image.open(f).convert("RGB") for f, _ in files]
w, h = imgs[0].size
pad, lab = 10, 46
canvas = Image.new("RGB", (w*5 + pad*6, h + lab + pad*2), (255, 255, 255))
d = ImageDraw.Draw(canvas)
try:
    fnt = ImageFont.truetype("arialbd.ttf", 34)
except Exception:
    fnt = ImageFont.load_default()
for i, (im, (_, t)) in enumerate(zip(imgs, files)):
    x = pad + i*(w+pad)
    canvas.paste(im, (x, lab))
    d.text((x + 12, 10), t, fill=(40, 40, 40), font=fnt)
canvas.save("fig6_snapshots.png", dpi=(300, 300))
print("saved fig5_contacts.png and fig6_snapshots.png", canvas.size)
