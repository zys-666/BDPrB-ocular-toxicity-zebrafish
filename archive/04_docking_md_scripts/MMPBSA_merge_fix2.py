# -*- coding: utf-8 -*-
"""Merge the PBC-corrected frames into the MM-PBSA statistics (robust version)."""
import glob, os, re, statistics as st
import numpy as np

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA"
num = re.compile(r"[-+]?\d+\.\d+")
bad_t = [58.4, 58.7, 59.3, 60.5, 60.8, 87.0, 87.6, 87.9, 88.2, 89.7, 95.3, 96.5]


def parse_total(path):
    out = []
    for line in open(path):
        if not line.startswith("_fix"):
            continue
        label = line.split()[0]
        t = float(label.split("~")[1].replace("ns", ""))
        rest = line[len(label):]
        p = rest.split("|")
        b = num.findall(p[0]); mm = num.findall(p[1]); c = num.findall(p[2])
        out.append(dict(t=t, bind=float(b[0]), bind_dh=float(b[1]), MM=float(mm[0]),
                        PB=float(mm[2]), SA=float(mm[3]), COU=float(c[0]), VDW=float(c[2])))
    return out

fix = []
for f in glob.glob(os.path.join(base, "fix12", "fix*", "_fix*~MMPBSA.dat")):
    fix += parse_total(f)
# dedupe by time, keep the run where the ligand is actually bound (MM < -10)
fix = [r for r in fix if r["MM"] < -10]
uniq = {}
for r in fix:
    k = round(r["t"], 1)
    if k not in uniq or r["MM"] < uniq[k]["MM"]:
        uniq[k] = r
fix = [uniq[k] for k in sorted(uniq)]
print("unique fixed frames:", len(fix))

old = np.loadtxt(os.path.join(base, "_pid~MMPBSA.dat"), comments="#")
rows = []
for i, t in enumerate(old[:, 0]):
    if any(abs(t - b) < 0.05 for b in bad_t):
        continue
    rows.append(dict(t=t, bind=old[i, 1], bind_dh=old[i, 2], MM=old[i, 3],
                     COU=old[i, 4], VDW=old[i, 5], PB=old[i, 6], SA=old[i, 7]))
rows += fix
rows.sort(key=lambda r: r["t"])
print("combined frames:", len(rows))


def rep(v):
    m = st.mean(v); sd = st.stdev(v); sem = sd / len(v) ** 0.5
    return m, sem, sd


print("\n=== MM-PBSA after PBC fix (%d frames) ===" % len(rows))
for key, name in [("bind_dh", "dG(PB, DH)"), ("MM", "dE_MM"), ("COU", "dE_COU"),
                  ("VDW", "dE_VDW"), ("PB", "dG_PB"), ("SA", "dG_SA")]:
    m, sem, sd = rep([r[key] for r in rows])
    print("  %-12s %9.2f +/- %5.2f (SD %5.2f) kJ/mol = %8.2f +/- %4.2f kcal/mol" % (
        name, m, sem, sd, m / 4.184, sem / 4.184))

out = os.path.join(base, "_pid~MMPBSA_fixed.dat")
with open(out, "w") as f:
    f.write("#Frame  dG_bind(kJ/mol)  dG_bind_DH(kJ/mol)  dE_MM  dE_COU  dE_VDW  dG_PB  dG_SA\n")
    for r in rows:
        f.write("%8.1f  %10.3f  %10.3f  %9.3f  %8.3f  %9.3f  %8.3f  %7.3f\n" % (
            r["t"], r["bind"], r["bind_dh"], r["MM"], r["COU"], r["VDW"], r["PB"], r["SA"]))
print("wrote", out)


def parse_res(path):
    out = []
    for line in open(path):
        if not (line.startswith("_seg") or line.startswith("_fix")):
            continue
        label = line.split()[0]
        t = float(label.split("~")[1].replace("ns", ""))
        vals = num.findall(line[len(label):])
        out.append((t, [float(vals[2 * i + 1]) for i in range(len(vals) // 2)]))
    return out

frames = []
for f in glob.glob(os.path.join(base, "par*", "seg*", "_seg*~res_MMPBSA.dat")):
    frames += parse_res(f)
print("\nseg frames parsed:", len(frames))
fixres = []
for f in glob.glob(os.path.join(base, "fix12", "fix*", "_fix*~res_MMPBSA.dat")):
    fixres += parse_res(f)
print("fix res frames parsed:", len(fixres))

keep = [fr for fr in frames if not any(abs(fr[0] - b) < 0.05 for b in bad_t)]
fixkeep = []
seen = set()
for fr in fixres:
    k = round(fr[0], 1)
    if k in seen or k not in {round(r["t"], 1) for r in fix}:
        continue
    seen.add(k)
    fixkeep.append(fr)
print("kept: %d bad-removed + %d fixed" % (len(keep), len(fixkeep)))
allf = keep + fixkeep
arr = np.array([fr[1] for fr in allf])
means = arr.mean(0)

order = []
for line in open(os.path.join(base, "_pid~res.dat")):
    if line.startswith("#"):
        continue
    p = line.split()
    if len(p) >= 2:
        order.append((p[0], p[1]))

with open(os.path.join(base, "_pid~res_fixed.dat"), "w") as f:
    f.write("#mol res  mean_DH(kJ/mol)  SD(kJ/mol)\n")
    for i, (mol, res) in enumerate(order):
        if i >= arr.shape[1]:
            break
        f.write("%-4s %-8s %10.3f %10.3f\n" % (mol, res, means[i], arr[:, i].std(ddof=1)))
print("wrote _pid~res_fixed.dat (%d frames, %d columns)" % (arr.shape[0], arr.shape[1]))
