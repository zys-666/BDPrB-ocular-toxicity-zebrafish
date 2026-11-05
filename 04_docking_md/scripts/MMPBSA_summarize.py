# -*- coding: utf-8 -*-
import glob, re, os, statistics as st

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA"
seg_files = glob.glob(os.path.join(base, "par*", "seg*", "_seg*~MMPBSA.dat"))
print("segment files:", len(seg_files))

num = re.compile(r"[-+]?\d+\.\d+")
rows = []
for f in seg_files:
    for line in open(f):
        if not line.startswith("_seg"):
            continue
        label = line.split()[0]
        t = float(label.split("~")[1].replace("ns", ""))
        parts = line.split("|")
        b = num.findall(parts[0])
        mm = num.findall(parts[1])
        cou = num.findall(parts[2])
        rows.append(dict(t=t, bind=float(b[0]), bind_dh=float(b[1]),
                         MM=float(mm[0]), MM_dh=float(mm[1]), PB=float(mm[2]), SA=float(mm[3]),
                         COU=float(cou[0]), COU_dh=float(cou[1]), VDW=float(cou[2])))
rows.sort(key=lambda r: r["t"])
print("frames:", len(rows), "time range:", rows[0]["t"], rows[-1]["t"])

def stats(vals):
    m = st.mean(vals); sd = st.stdev(vals); sem = sd/len(vals)**0.5
    return m, sd, sem

print()
print("=== MM-PBSA (last 50 ns, %d frames) ===" % len(rows))
for key, name in [("bind","dG_PB (kJ/mol)"), ("bind_dh","dG_PB+DH (kJ/mol)"),
                  ("MM","dE_MM"), ("COU","dE_COU"), ("VDW","dE_VDW"),
                  ("PB","dG_PB(solv)"), ("SA","dG_SA")]:
    vals = [r[key] for r in rows]
    m, sd, sem = stats(vals)
    print("  %-16s %9.2f  +/- %6.2f (SD %6.2f)  = %7.2f kcal/mol" % (name, m, sem, sd, m/4.184))

# write clean merged file
out = os.path.join(base, "_pid~MMPBSA.dat")
with open(out, "w") as f:
    f.write("#Frame  dG_bind(kJ/mol)  dG_bind_DH  dE_MM  dE_COU  dE_VDW  dG_PB  dG_SA\n")
    for r in rows:
        f.write("%8.1f  %10.3f  %10.3f  %9.3f  %8.3f  %9.3f  %8.3f  %7.3f\n" % (
            r["t"], r["bind"], r["bind_dh"], r["MM"], r["COU"], r["VDW"], r["PB"], r["SA"]))
print("wrote", out)

# per-residue
res = []
for line in open(os.path.join(base, "_pid~res.dat")):
    if line.startswith("#"): continue
    p = line.split()
    if len(p) < 3: continue
    label = p[0] + p[1]
    vals = num.findall(line)
    if len(vals) < 4: continue
    res.append((p[0], p[1], float(vals[0]), float(vals[1])))
prot = [x for x in res if x[0] != "Lig"]
prot.sort(key=lambda x: x[2])
print()
print("=== top 15 favourable residues (per-residue MM-PBSA, kJ/mol) ===")
for mol, r, v, vdh in prot[:15]:
    print("   %-8s %8.2f  (with DH %8.2f)" % (r, v, vdh))
print()
print("=== top 5 unfavourable ===")
for mol, r, v, vdh in prot[-5:]:
    print("   %-8s %8.2f" % (r, v))
