# -*- coding: utf-8 -*-
import glob, re, os, statistics as st

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA"
num = re.compile(r"[-+]?\d+\.\d+")
rows = []
for f in glob.glob(os.path.join(base, "par*", "seg*", "_seg*~MMPBSA.dat")):
    for line in open(f):
        if not line.startswith("_seg"):
            continue
        label = line.split()[0]
        t = float(label.split("~")[1].replace("ns", ""))
        rest = line[len(label):]
        parts = rest.split("|")
        b = num.findall(parts[0]); mm = num.findall(parts[1]); cou = num.findall(parts[2])
        rows.append(dict(t=t, bind=float(b[0]), bind_dh=float(b[1]),
                         MM=float(mm[0]), PB=float(mm[2]), SA=float(mm[3]),
                         COU=float(cou[0]), COU_dh=float(cou[1]), VDW=float(cou[2])))
rows.sort(key=lambda r: r["t"])
print("frames:", len(rows), "time:", rows[0]["t"], "-", rows[-1]["t"])

def stats(v):
    return st.mean(v), st.stdev(v)/len(v)**0.5, st.stdev(v)

print()
print("=== MM-PBSA, last 50 ns, %d frames ===" % len(rows))
for key, name in [("bind_dh","dG(PB, with DH)"), ("bind","dG(PB, no DH)"),
                  ("MM","dE_MM"), ("COU","dE_COU"), ("VDW","dE_VDW"),
                  ("PB","dG_PB(solv)"), ("SA","dG_SA")]:
    m, sem, sd = stats([r[key] for r in rows])
    print("  %-18s %9.2f +/- %5.2f (SD %5.2f) kJ/mol = %8.2f +/- %4.2f kcal/mol" % (name, m, sem, sd, m/4.184, sem/4.184))

out = os.path.join(base, "_pid~MMPBSA.dat")
with open(out, "w") as f:
    f.write("#Frame  dG_bind(kJ/mol)  dG_bind_DH(kJ/mol)  dE_MM  dE_COU  dE_VDW  dG_PB  dG_SA\n")
    for r in rows:
        f.write("%8.1f  %10.3f  %10.3f  %9.3f  %8.3f  %9.3f  %8.3f  %7.3f\n" % (
            r["t"], r["bind"], r["bind_dh"], r["MM"], r["COU"], r["VDW"], r["PB"], r["SA"]))
print("\nwrote", out)

# sanity: component sum
m_bind = st.mean([r["bind_dh"] for r in rows]); m_mm = st.mean([r["MM"] for r in rows])
m_pb = st.mean([r["PB"] for r in rows]); m_sa = st.mean([r["SA"] for r in rows])
print("check MM+PB+SA = %.2f vs dG(DH) = %.2f" % (m_mm+m_pb+m_sa, m_bind))
