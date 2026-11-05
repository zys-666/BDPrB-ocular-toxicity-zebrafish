# -*- coding: utf-8 -*-
"""Summarise pocket RMSF and pocket SASA for the three systems."""
import numpy as np, os

root = r"D:\user data\Desktop\斑马鱼\5.dockMD"
dirs = {
    "apo":   os.path.join(root, "controls", "apo"),
    "BDPrB": os.path.join(root, "MD", "03_production"),
    "KyA33": os.path.join(root, "controls", "kya33"),
}
pocket = [138, 141, 142, 145, 146, 188, 191, 192, 195, 262, 314, 319, 320, 321,
          475, 477, 478, 479, 483]

print("%-8s %14s %14s" % ("system", "pocket RMSF(A)", "pocket SASA(nm2)"))
for name, d in dirs.items():
    r = np.loadtxt(os.path.join(d, "rmsf_all.xvg"), comments=("#", "@"))
    # first column is residue number, second RMSF in nm
    res = r[:, 0].astype(int)
    rmsf = r[:, 1] * 10
    vals = [rmsf[res == x][0] for x in pocket if (res == x).any()]
    s = np.loadtxt(os.path.join(d, "pocket_sasa.xvg"), comments=("#", "@"))
    sasa = s[:, 1].mean()
    print("%-8s %14.2f %14.2f" % (name, np.mean(vals), sasa))

print("\nper-residue RMSF (A):")
print("%-8s %8s %8s %8s" % ("resid", "apo", "BDPrB", "KyA33"))
data = {}
for name, d in dirs.items():
    r = np.loadtxt(os.path.join(d, "rmsf_all.xvg"), comments=("#", "@"))
    res = r[:, 0].astype(int)
    data[name] = (res, r[:, 1] * 10)
for x in pocket:
    row = []
    for name in ("apo", "BDPrB", "KyA33"):
        res, v = data[name]
        row.append(v[res == x][0] if (res == x).any() else float("nan"))
    print("%-8d %8.2f %8.2f %8.2f" % (x, *row))
