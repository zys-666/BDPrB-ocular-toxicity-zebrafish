# -*- coding: utf-8 -*-
import re, os
base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA"
num = re.compile(r"[-+]?\d+\.\d+")
rows = []
for line in open(os.path.join(base, "_pid~res.dat")):
    if line.startswith("#"): continue
    p = line.split()
    if len(p) < 3: continue
    v = num.findall(line)
    rows.append((p[0], p[1], float(v[0]), float(v[1])))
d = {r[1]: r[2] for r in rows}
print("ligand row:", d.get("1LIG"))
pocket = ["138VAL","188PHE","192MET","195TRP","314PHE","319CYS","320CYS","321THR","475ASN","477MET","478SER","479CYS","483PHE","145LYS","191VAL","262THR"]
print()
print("口袋残基贡献 (kJ/mol):")
for k in pocket:
    if k in d: print("   %-8s %7.2f  (%6.2f kcal/mol)" % (k, d[k], d[k]/4.184))
prot_sum = sum(r[2] for r in rows if r[0] == "Pro")
print("\nprotein contributions sum: %.2f kJ/mol" % prot_sum)
