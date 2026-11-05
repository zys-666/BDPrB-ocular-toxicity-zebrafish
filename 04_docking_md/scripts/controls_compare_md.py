# -*- coding: utf-8 -*-
"""Compare protein / ligand / pocket behaviour across apo, BDPrB-holo and KyA33-holo."""
import numpy as np, os

root = r"D:\user data\Desktop\斑马鱼\5.dockMD"
systems = {
    "apo":   os.path.join(root, "controls", "apo"),
    "BDPrB": os.path.join(root, "MD", "04_analysis"),
    "KyA33": os.path.join(root, "controls", "kya33"),
}

def load(p):
    d = np.loadtxt(p, comments=("#", "@"))
    t = d[:, 0]
    return t, d

print("%-8s %10s %10s %10s %10s" % ("system", "RMSD(A)", "last50(A)", "Rg(A)", "ligRMSD(A)"))
for name, d in systems.items():
    t, r = load(os.path.join(d, "rmsd_bb.xvg"))
    rms = r[:, 1] * 10
    _, g = load(os.path.join(d, "gyrate.xvg"))
    rg = g[:, 1] * 10
    m = t >= 50
    lig = ""
    lp = os.path.join(d, "rmsd_lig.xvg")
    if name != "apo" and os.path.exists(lp):
        _, l = load(lp)
        lig = "%10.2f" % (l[m, 1].mean() * 10)
    print("%-8s %10.2f %10.2f %10.2f %s" % (name, rms.mean(), rms[m].mean(), rg.mean(), lig))
