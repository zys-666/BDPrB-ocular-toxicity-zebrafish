# -*- coding: utf-8 -*-
"""Per-residue ligand-contact occupancy from a multi-model complex PDB."""
import numpy as np, sys, os

CUT = 4.0

def analyse(path, out_csv):
    prot_res = []      # residue ids in atom order
    frames = []
    cur_p, cur_l = [], []
    cur_ids = []
    model = False
    for line in open(path):
        if line.startswith("MODEL"):
            model = True; cur_p, cur_l, cur_ids = [], [], []
        elif line.startswith("ENDMDL"):
            frames.append((np.array(cur_p), np.array(cur_l), cur_ids))
            model = False
        elif line.startswith(("ATOM", "HETATM")):
            ch = line[21]
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if ch == "A" or ch == " ":
                cur_p.append(xyz); cur_ids.append(int(line[22:26]))
            else:
                cur_l.append(xyz)
    if not frames:
        print("no MODEL records:", path); return None
    print("%s: %d frames, %d protein atoms, %d ligand atoms" % (
        os.path.basename(path), len(frames), len(frames[0][0]), len(frames[0][1])))

    ids = np.array(frames[0][2])
    uniq = sorted(set(ids.tolist()))
    occ = {r: 0 for r in uniq}
    for P, L, _ in frames:
        if len(L) == 0:
            continue
        d = np.sqrt(((P[:, None, :] - L[None, :, :]) ** 2).sum(-1))   # (natom, nlig)
        dmin = d.min(1)
        # per residue min
        for r in uniq:
            sel = ids == r
            if dmin[sel].min() <= CUT:
                occ[r] += 1
    n = len(frames)
    with open(out_csv, "w") as f:
        f.write("#resid\toccupancy_pct\n")
        for r in uniq:
            f.write("%d\t%.1f\n" % (r, 100.0 * occ[r] / n))
    print("  wrote", out_csv)
    return {r: 100.0 * occ[r] / n for r in uniq}

a = analyse(r"D:\user data\Desktop\斑马鱼\5.dockMD\MD\03_production\traj167.pdb",
            r"D:\user data\Desktop\斑马鱼\5.dockMD\MD\03_production\contacts_BDPrB.csv")
b = analyse(r"D:\user data\Desktop\斑马鱼\5.dockMD\controls\kya33\traj167.pdb",
            r"D:\user data\Desktop\斑马鱼\5.dockMD\controls\kya33\contacts_KyA33.csv")
