# -*- coding: utf-8 -*-
"""Per-residue ligand-contact occupancy from a multi-model complex PDB."""
import numpy as np, os

CUT = 4.0

def analyse(path, out_csv):
    frames, cur_p, cur_l, cur_ids = [], [], [], []
    for line in open(path):
        if line.startswith("MODEL"):
            cur_p, cur_l, cur_ids = [], [], []
        elif line.startswith("ENDMDL"):
            frames.append((np.array(cur_p), np.array(cur_l), list(cur_ids)))
        elif line.startswith(("ATOM", "HETATM")):
            ch = line[21]
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if ch == "A" or ch == " ":
                cur_p.append(xyz); cur_ids.append(int(line[22:26]))
            else:
                cur_l.append(xyz)
    if not frames:
        print("no frames:", path); return None
    print("%s: %d frames, %d protein atoms, %d ligand atoms" % (
        os.path.basename(path), len(frames), len(frames[0][0]), len(frames[0][1])))

    ids = np.array(frames[0][2])
    uniq = sorted(set(ids.tolist()))
    idx = {r: np.where(ids == r)[0] for r in uniq}
    occ = np.zeros(len(uniq))
    for P, L, _ in frames:
        if len(L) == 0:
            continue
        d = np.sqrt(((P[:, None, :] - L[None, :, :]) ** 2).sum(-1))
        dmin = d.min(1)
        for k, r in enumerate(uniq):
            if dmin[idx[r]].min() <= CUT:
                occ[k] += 1
    n = len(frames)
    with open(out_csv, "w") as f:
        f.write("#resid\toccupancy_pct\n")
        for k, r in enumerate(uniq):
            f.write("%d\t%.1f\n" % (r, 100.0 * occ[k] / n))
    print("  wrote", out_csv)
    return dict(zip(uniq, 100.0 * occ / n))

a = analyse("../MD/03_production/traj167.pdb", "../MD/03_production/contacts_BDPrB.csv")
b = analyse("../controls/kya33/traj167.pdb", "../controls/kya33/contacts_KyA33.csv")
if a and b:
    top_a = sorted(a.items(), key=lambda kv: -kv[1])[:12]
    top_b = sorted(b.items(), key=lambda kv: -kv[1])[:12]
    print("BDPrB top:", [(r, round(v)) for r, v in top_a])
    print("KyA33 top:", [(r, round(v)) for r, v in top_b])
