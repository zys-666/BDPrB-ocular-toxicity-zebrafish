# -*- coding: utf-8 -*-
import numpy as np, re, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA"
d = np.loadtxt(os.path.join(base, "_pid~MMPBSA.dat"), comments="#")
t, dg = d[:,0], d[:,1]

num = re.compile(r"[-+]?\d+\.\d+")
res = []
for line in open(os.path.join(base, "_pid~res.dat")):
    if line.startswith("#"): continue
    p = line.split()
    if len(p) < 3 or p[0] != "Pro": continue
    v = num.findall(line)
    res.append((p[1], float(v[0])))
res.sort(key=lambda x: x[1])

fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
ax[0].plot(t, dg, lw=0.8, color="tab:blue")
ax[0].axhline(dg.mean(), color="tab:red", ls="--", lw=1, label="mean %.1f kJ/mol" % dg.mean())
ax[0].set_xlabel("Time (ns)"); ax[0].set_ylabel("dG_bind (kJ/mol)")
ax[0].set_title("MM-PBSA binding free energy")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)

top = res[:12][::-1]
ax[1].barh([x[0] for x in top], [x[1] for x in top], color="tab:green")
ax[1].set_xlabel("Contribution (kJ/mol)")
ax[1].set_title("Top favourable residues")
ax[1].grid(alpha=0.3, axis="x")
plt.tight_layout()
out = os.path.join(base, "mmpbsa_summary.png")
plt.savefig(out, dpi=300)
print("saved", out)
print("dG mean %.1f, SD %.1f, first-half %.1f, second-half %.1f" % (
    dg.mean(), dg.std(), dg[:len(dg)//2].mean(), dg[len(dg)//2:].mean()))
