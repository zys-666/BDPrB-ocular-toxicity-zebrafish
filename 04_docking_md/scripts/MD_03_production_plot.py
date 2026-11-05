import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def load(p):
    return np.loadtxt(p, comments=("#", "@"))

rmsd = load(r"..\04_analysis\rmsd_bb.xvg")
lig = load(r"..\04_analysis\rmsd_lig.xvg")
rg = load(r"..\04_analysis\gyrate.xvg")
c320 = load(r"c320_lig.xvg")
mind = load(r"mindist.xvg")

fig, ax = plt.subplots(2, 2, figsize=(10, 7))
ax[0,0].plot(rmsd[:,0], rmsd[:,1]*10, lw=1, color="tab:blue")
ax[0,0].set_xlabel("Time (ns)"); ax[0,0].set_ylabel("Backbone RMSD (A)")
ax[0,0].set_title("Protein backbone RMSD")
ax[0,0].set_ylim(0, 8)

ax[0,1].plot(lig[:,0], lig[:,1]*10, lw=1, color="tab:red")
ax[0,1].set_xlabel("Time (ns)"); ax[0,1].set_ylabel("Ligand RMSD (A)")
ax[0,1].set_title("BDPrB RMSD (fit on protein)")
ax[0,1].set_ylim(0, 8)

ax[1,0].plot(rg[:,0]/1000.0, rg[:,1]*10, lw=1, color="tab:green")
ax[1,0].set_xlabel("Time (ns)"); ax[1,0].set_ylabel("Rg (A)")
ax[1,0].set_title("Radius of gyration")

d = c320[:,1:].min(axis=1)*10
ax[1,1].plot(c320[:,0]/1000.0, d, lw=0.8, color="tab:purple", label="min dist to Cys320 SG")
ax[1,1].plot(mind[:,0]/1000.0, mind[:,1]*10, lw=0.8, color="tab:gray", alpha=0.7, label="min dist to protein")
ax[1,1].set_xlabel("Time (ns)"); ax[1,1].set_ylabel("Distance (A)")
ax[1,1].set_title("Ligand-pocket distances")
ax[1,1].legend(fontsize=8)
for a in ax.flat:
    a.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(r"..\04_analysis\md_summary.png", dpi=300)
print("saved md_summary.png")

