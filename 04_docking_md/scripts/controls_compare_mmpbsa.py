# -*- coding: utf-8 -*-
"""Final comparison: BDPrB vs KyA33 (MM-PBSA + docking + MD)."""
import numpy as np, os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

root = r"D:\user data\Desktop\斑马鱼\5.dockMD"
files = {
    "BDPrB": os.path.join(root, "MMPBSA", "_pid~MMPBSA_fixed.dat"),
    "KyA33": os.path.join(root, "MMPBSA_kya33", "_pid~MMPBSA.dat"),
}
data = {}
for k, f in files.items():
    d = np.loadtxt(f, comments="#")
    data[k] = d
    print("%-6s %4d frames  dG = %.1f +/- %.1f (SD %.1f) kJ/mol = %.2f kcal/mol" % (
        k, len(d), d[:, 2].mean(), d[:, 2].std(ddof=1) / len(d) ** 0.5, d[:, 2].std(ddof=1),
        d[:, 2].mean() / 4.184))

print()
print("%-12s %10s %10s" % ("component", "BDPrB", "KyA33"))
names = ["dE_MM", "dE_coul", "dE_vdW", "dG_PB", "dG_SA", "dG_bind"]
cols = [3, 4, 5, 6, 7, 2]
for nm, c in zip(names, cols):
    a = data["BDPrB"][:, c].mean(); b = data["KyA33"][:, c].mean()
    print("%-12s %10.1f %10.1f" % (nm, a, b))

fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for k, col in (("BDPrB", "tab:blue"), ("KyA33", "tab:orange")):
    ax[0].hist(data[k][:, 2], bins=25, alpha=0.55, label="%s (%.0f kJ/mol)" % (k, data[k][:, 2].mean()),
               color=col)
ax[0].set_xlabel("$\\Delta G_{bind}$ (kJ/mol)")
ax[0].set_ylabel("count")
ax[0].set_title("MM-PBSA distributions")
ax[0].legend(fontsize=8)

x = np.arange(len(names) - 1)
w = 0.38
ax[1].bar(x - w/2, [data["BDPrB"][:, c].mean() for c in cols[:-1]], w, label="BDPrB", color="tab:blue")
ax[1].bar(x + w/2, [data["KyA33"][:, c].mean() for c in cols[:-1]], w, label="KyA33", color="tab:orange")
ax[1].set_xticks(x); ax[1].set_xticklabels(names[:-1], rotation=20)
ax[1].axhline(0, color="k", lw=0.6)
ax[1].set_ylabel("kJ/mol")
ax[1].set_title("Energy decomposition")
ax[1].legend(fontsize=8)
plt.tight_layout()
out = os.path.join(root, "controls", "mmpbsa_comparison.png")
plt.savefig(out, dpi=300)
print("\nsaved", out)
