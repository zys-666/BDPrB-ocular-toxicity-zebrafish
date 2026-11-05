# -*- coding: utf-8 -*-
import numpy as np, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA"
d = np.loadtxt(os.path.join(base, "_pid~MMPBSA_fixed.dat"), comments="#")
t, dg = d[:, 0], d[:, 2]          # column 2 = dG with Debye-Huckel correction

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ax[0].plot(t, dg, lw=0.8, color="tab:blue")
ax[0].axhline(dg.mean(), color="tab:red", ls="--", lw=1,
              label="mean %.1f $\\pm$ %.1f kJ/mol" % (dg.mean(), dg.std(ddof=1) / len(dg) ** 0.5))
ax[0].set_xlabel("Time (ns)")
ax[0].set_ylabel("$\\Delta G_{bind}$ (kJ/mol)")
ax[0].set_title("MM-PBSA binding free energy (168 frames, with DH)")
ax[0].legend(fontsize=8)
ax[0].grid(alpha=0.3)

res = []
for line in open(os.path.join(base, "_pid~res_fixed.dat")):
    if line.startswith("#"):
        continue
    p = line.split()
    res.append((p[0], p[1], float(p[2]), float(p[3])))
prot = [r for r in res if r[0] == "Pro"]
prot.sort(key=lambda x: x[2])
top = prot[:12][::-1]
ax[1].barh([x[1] for x in top], [x[2] for x in top], color="tab:green")
ax[1].set_xlabel("Contribution (kJ/mol)")
ax[1].set_title("Top favourable residues")
ax[1].grid(alpha=0.3, axis="x")
plt.tight_layout()
plt.savefig(os.path.join(base, "mmpbsa_summary.png"), dpi=300)
print("mean %.2f  sem %.2f  sd %.2f" % (dg.mean(), dg.std(ddof=1) / len(dg) ** 0.5, dg.std(ddof=1)))
