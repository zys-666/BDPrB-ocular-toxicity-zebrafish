# -*- coding: utf-8 -*-
"""Export the docking best poses as clean PDB / SDF (PDBQT keeps no bond orders)."""
import os, subprocess, shutil

root = r"D:\user data\Desktop\斑马鱼\5.dockMD"
out = os.path.join(root, "dock", "docking_results")
os.makedirs(out, exist_ok=True)
obabel = r"D:\python313\Scripts\obabel.exe"

# source mol2 files that already carry the docked coordinates + correct bond orders
jobs = [
    ("BDPrB_pose1", os.path.join(root, "MD", "00_prep", "lig_mapped.mol2")),
    ("KyA33_zebrafish_pose1", os.path.join(root, "controls", "kya33", "kya33_pose.mol2")),
    ("KyA33_redock_pose1", os.path.join(root, "dock", "recontrol", "pose1.mol2")),
]
for name, src in jobs:
    if not os.path.exists(src):
        print("MISSING", src); continue
    for ext in ("pdb", "sdf"):
        dst = os.path.join(out, "%s.%s" % (name, ext))
        subprocess.run([obabel, src, "-O", dst], capture_output=True)
    print("wrote", name, ".pdb / .sdf")

# protein + ligand complexes for PyMOL viewing
def merge(prot, lig, dst, chain="L"):
    lines = [l for l in open(prot) if l.startswith(("ATOM", "HETATM"))]
    with open(dst, "w") as f:
        for l in lines:
            f.write(l if l.endswith("\n") else l + "\n")
        f.write("TER\n")
        n = len(lines)
        for l in open(lig):
            if l.startswith(("ATOM", "HETATM")):
                n += 1
                f.write(l[:21] + chain + l[22:])
        f.write("END\n")
    print("wrote", dst)

merge(os.path.join(root, "dock", "zebrafish_ALDH1A2_trim26-518.pdb"),
      os.path.join(out, "BDPrB_pose1.pdb"),
      os.path.join(out, "complex_BDPrB.pdb"))
merge(os.path.join(root, "dock", "zebrafish_ALDH1A2_trim26-518.pdb"),
      os.path.join(out, "KyA33_zebrafish_pose1.pdb"),
      os.path.join(out, "complex_KyA33.pdb"))
