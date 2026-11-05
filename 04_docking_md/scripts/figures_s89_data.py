# -*- coding: utf-8 -*-
"""Prepare data for Figures S8 / S9."""
import numpy as np, os, json

BASE = r"D:\user data\Desktop\斑马鱼"
OUT = r"D:\user data\Desktop\斑马鱼\5.dockMD\figures"

# ---------- S8: pLDDT from the full AF2 model ----------
p = os.path.join(BASE, "14.补充材料", "AF-Q90XS8-F1-model_v6.pdb")
res, plddt = [], []
for line in open(p):
    if line.startswith("ATOM") and line[12:16].strip() == "CA":
        res.append(int(line[22:26]))
        plddt.append(float(line[60:66]))
res = np.array(res); plddt = np.array(plddt)
np.savez(os.path.join(OUT, "s8_plddt.npz"), res=res, plddt=plddt)
print("pLDDT: %d residues, range %d-%d, mean %.1f" % (len(res), res.min(), res.max(), plddt.mean()))
print("  1-25 mean %.1f (min %.1f max %.1f)" % (plddt[res <= 25].mean(), plddt[res <= 25].min(), plddt[res <= 25].max()))
print("  26-518 mean %.1f" % plddt[(res >= 26)].mean())
print("  pLDDT>70: %.1f%%" % (100.0 * (plddt >= 70).mean()))

# ---------- S9: contact heatmap (residue x time) ----------
def heatmap(path, out_npz, cut=4.0):
    frames, cur_p, cur_l, cur_ids = [], [], [], []
    for line in open(path):
        if line.startswith("MODEL"):
            cur_p, cur_l, cur_ids = [], [], []
        elif line.startswith("ENDMDL"):
            frames.append((np.array(cur_p), np.array(cur_l), np.array(cur_ids)))
        elif line.startswith(("ATOM", "HETATM")):
            ch = line[21]
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            if ch in ("A", " "):
                cur_p.append(xyz); cur_ids.append(int(line[22:26]))
            else:
                cur_l.append(xyz)
    ids = frames[0][2]
    uniq = np.unique(ids)
    mat = np.zeros((len(uniq), len(frames)), dtype=np.float32)
    for fi, (P, L, _) in enumerate(frames):
        d = np.linalg.norm(P[:, None, :] - L[None, :, :], axis=-1).min(1)
        for k, r in enumerate(uniq):
            mat[k, fi] = 1.0 if d[ids == r].min() <= cut else 0.0
    np.savez(out_npz, res=uniq, mat=mat)
    print("%s: %d residues x %d frames, mean contacts/frame %.1f" % (
        os.path.basename(out_npz), mat.shape[0], mat.shape[1], mat.sum(0).mean()))

heatmap(os.path.join(BASE, "5.dockMD", "MD", "03_production", "traj167.pdb"),
        os.path.join(OUT, "s9_heat_BDPrB.npz"))
heatmap(os.path.join(BASE, "5.dockMD", "controls", "kya33", "traj167.pdb"),
        os.path.join(OUT, "s9_heat_KyA33.npz"))
print("done")
