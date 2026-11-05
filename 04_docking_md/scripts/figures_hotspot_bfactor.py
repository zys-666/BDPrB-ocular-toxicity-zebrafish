# -*- coding: utf-8 -*-
"""Write complex PDBs whose B-factor = per-residue MM-PBSA contribution (for hot-spot surfaces)."""
import re, os

def read_res(path, dh=False):
    d = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.split()
        if len(p) < 3 or p[0] == "Lig":
            continue
        nums = re.findall(r"[-+]?\d+\.\d+", line)
        if not nums:
            continue
        d[int(p[1][:-3])] = float(nums[0])
    return d

def write(src, dst, vals):
    n = 0
    with open(dst, "w") as f:
        for l in open(src):
            if not l.startswith(("ATOM", "HETATM")):
                continue
            n += 1
            rid = int(l[22:26])
            b = vals.get(rid, 0.0)
            f.write(l[:60] + "%6.2f" % b + l[66:].rstrip() + "\n")
    print("wrote", dst, n, "atoms, residues with data:", len(vals))

write("../dock/docking_results/complex_BDPrB.pdb",
      "hotspot_BDPrB_bfactor.pdb",
      read_res("../MMPBSA/_pid~res_fixed.dat"))
write("../dock/docking_results/complex_KyA33.pdb",
      "hotspot_KyA33_bfactor.pdb",
      read_res("../MMPBSA_kya33/_pid~res.dat"))
