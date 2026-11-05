# -*- coding: utf-8 -*-
"""Aggregate the KyA33 MM-PBSA segment results (the driver's per-frame merge produced 0 rows)."""
import glob, os, re, statistics as st
import numpy as np

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA_kya33"
num = re.compile(r"[-+]?\d+\.\d+")

rows = []
for f in glob.glob(os.path.join(base, "par*", "seg*", "_seg*~MMPBSA.dat")):
    for line in open(f):
        if not line.startswith("_seg"):
            continue
        label = line.split()[0]
        t = float(label.split("~")[1].replace("ns", ""))
        rest = line[len(label):]
        p = rest.split("|")
        b = num.findall(p[0]); mm = num.findall(p[1]); c = num.findall(p[2])
        rows.append(dict(t=t, bind=float(b[0]), bind_dh=float(b[1]), MM=float(mm[0]),
                         PB=float(mm[2]), SA=float(mm[3]), COU=float(c[0]), VDW=float(c[2])))
rows.sort(key=lambda r: r["t"])
print("frames:", len(rows), "time %.1f - %.1f ns" % (rows[0]["t"], rows[-1]["t"]))

bad = [r for r in rows if r["MM"] > -10]
print("artifact frames (MM ~ 0):", len(bad), [r["t"] for r in bad])
good = [r for r in rows if r["MM"] <= -10]


def rep(v):
    m = st.mean(v); sd = st.stdev(v); sem = sd / len(v) ** 0.5
    return m, sem, sd


print("\n=== KyA33 MM-PBSA (%d frames) ===" % len(good))
for key, name in [("bind_dh", "dG(PB, DH)"), ("MM", "dE_MM"), ("COU", "dE_COU"),
                  ("VDW", "dE_VDW"), ("PB", "dG_PB"), ("SA", "dG_SA")]:
    m, sem, sd = rep([r[key] for r in good])
    print("  %-12s %9.2f +/- %5.2f (SD %5.2f) kJ/mol = %8.2f +/- %4.2f kcal/mol" % (
        name, m, sem, sd, m / 4.184, sem / 4.184))

out = os.path.join(base, "_pid~MMPBSA.dat")
with open(out, "w") as f:
    f.write("#Frame  dG_bind(kJ/mol)  dG_bind_DH(kJ/mol)  dE_MM  dE_COU  dE_VDW  dG_PB  dG_SA\n")
    for r in good:
        f.write("%8.1f  %10.3f  %10.3f  %9.3f  %8.3f  %9.3f  %8.3f  %7.3f\n" % (
            r["t"], r["bind"], r["bind_dh"], r["MM"], r["COU"], r["VDW"], r["PB"], r["SA"]))
print("wrote", out)
