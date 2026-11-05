# -*- coding: utf-8 -*-
"""Protein (9R6N chain A) + KyA33 crystal pose (chain C) + redocked pose (chain D)."""
import os
import numpy as np

root = r"D:\user data\Desktop\斑马鱼\5.dockMD"
prot_p = os.path.join(root, "dock", "recontrol", "protein.pdb")
cry_p = os.path.join(root, "dock", "recontrol", "lig_ref.pdb")
red_p = os.path.join(root, "dock", "docking_results", "KyA33_redock_pose1.pdb")
out_p = os.path.join(root, "dock", "docking_results", "KyA33_crystal_vs_redock_in_protein.pdb")

ELS = {"C": "C", "N": "N", "O": "O", "S": "S", "H": "H", "F": "F"}
n = 0
with open(out_p, "w") as f:
    f.write("REMARK   KyA33 crystal pose (chain C) vs redocked pose (chain D) in ALDH1A2 (chain A)\n")
    f.write("REMARK   redock heavy-atom RMSD = 0.52 A (symmetry corrected, 23 atoms)\n")
    for path, chain in ((prot_p, "A"), (cry_p, "C"), (red_p, "D")):
        for l in open(path):
            if not l.startswith(("ATOM", "HETATM")):
                continue
            l = l.rstrip("\n").ljust(80)
            n += 1
            name = l[12:16].strip()
            el = l[76:78].strip()
            if not el:
                el = name[0] if name[0].isalpha() else "C"
            el = ELS.get(el.upper(), el.upper())
            f.write("%-6s%5d %-4s%1s%3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n" % (
                ("HETATM" if chain in ("C","D") else "ATOM"), n, name, " ",
                ("LIG" if chain in ("C","D") else l[17:20].strip()), chain,
                int(l[22:26]), float(l[30:38]), float(l[38:46]), float(l[46:54]), el))
        f.write("TER\n")
    f.write("END\n")
print("wrote", out_p, n, "atoms")

# quick sanity check
chains = {}
for l in open(out_p):
    if l.startswith(("ATOM", "HETATM")):
        chains[l[21]] = chains.get(l[21], 0) + 1
print("chains:", chains)
