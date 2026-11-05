# -*- coding: utf-8 -*-
"""Manual ligand-protein H-bond detection from multi-frame complex PDBs.

Criterion: H...A <= 2.5 A and angle(D-H...A) >= 120 deg.
Only inter-group (ligand <-> protein) bonds are counted.
"""
import numpy as np, sys, collections

HB_DIST = 2.5
HB_ANG = 120.0

ELEMS = {"C", "N", "O", "S", "F", "H", "P"}


def elem(name):
    n = name.strip()
    while n and n[-1].isdigit():
        n = n[:-1]
    if not n:
        return ""
    if len(n) >= 2 and n[:2].upper() in ("CL", "NA", "MG", "ZN", "FE"):
        return n[:2].upper()
    return n[0].upper()


def parse(path):
    frames = []
    cur = {"lig": [], "prot": []}
    for line in open(path):
        if line.startswith("MODEL"):
            cur = {"lig": [], "prot": []}
        elif line.startswith("ENDMDL"):
            frames.append(cur)
        elif line.startswith(("ATOM", "HETATM")):
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            rec = (line[12:16].strip(), line[17:20].strip(), line[22:26].strip(), xyz)
            (cur["lig"] if line[21] not in ("A", " ") else cur["prot"]).append(rec)
    if not frames:
        print("no MODEL records")
    return frames


def build_donors(atoms):
    """Return list of (D_index, H_index) plus arrays xyz/elem."""
    xyz = np.array([a[3] for a in atoms])
    el = [elem(a[0]) for a in atoms]
    hyd = [i for i, e in enumerate(el) if e == "H"]
    heavy = [i for i, e in enumerate(el) if e in ("N", "O", "S")]
    pairs = []
    if hyd and heavy:
        H = xyz[hyd]
        D = xyz[heavy]
        d = np.linalg.norm(H[:, None, :] - D[None, :, :], axis=-1)
        for k, hi in enumerate(hyd):
            j = int(np.argmin(d[k]))
            if d[k, j] <= 1.35:
                pairs.append((heavy[j], hi))
    return xyz, el, pairs


def analyse(path, tag, out_csv):
    frames = parse(path)
    n = len(frames)
    print("%s: %d frames, lig %d atoms, prot %d atoms" % (
        tag, n, len(frames[0]["lig"]), len(frames[0]["prot"])))

    lig_names = [a[0] for a in frames[0]["lig"]]
    prot_meta = [(a[1], a[2], a[0]) for a in frames[0]["prot"]]

    counts = np.zeros(n)
    occ = collections.Counter()

    for fi, fr in enumerate(frames):
        lxyz, lel, ldon = build_donors(fr["lig"])
        pxyz, pel, pdon = build_donors(fr["prot"])
        lig_acc = [i for i, e in enumerate(lel) if e in ("N", "O", "S")]
        prot_acc = [i for i, e in enumerate(pel) if e in ("N", "O", "S")]

        hits = []
        # ligand donor -> protein acceptor
        for di, hi in ldon:
            for ai in prot_acc:
                v1 = lxyz[di] - lxyz[hi]
                v2 = pxyz[ai] - lxyz[hi]
                r = np.linalg.norm(v2)
                if r > HB_DIST or r < 1e-3:
                    continue
                cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * r + 1e-9)
                if cosang < np.cos(np.radians(HB_ANG)):
                    continue
                hits.append(("L%s" % lig_names[di], prot_meta[ai][0] + prot_meta[ai][1]))
        # protein donor -> ligand acceptor
        for di, hi in pdon:
            for ai in lig_acc:
                v1 = pxyz[di] - pxyz[hi]
                v2 = lxyz[ai] - pxyz[hi]
                r = np.linalg.norm(v2)
                if r > HB_DIST or r < 1e-3:
                    continue
                cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * r + 1e-9)
                if cosang < np.cos(np.radians(HB_ANG)):
                    continue
                hits.append((prot_meta[di][0] + prot_meta[di][1], "L%s" % lig_names[ai]))
        counts[fi] = len(set(hits))
        for h in set(hits):
            occ[h] += 1

    with open(out_csv, "w") as f:
        f.write("#pair\toccupancy_pct\n")
        for pair, c in occ.most_common():
            f.write("%s-%s\t%.1f\n" % (pair[0], pair[1], 100.0 * c / n))
    print("  frames with >=1 H-bond: %d/%d" % (int((counts > 0).sum()), n))
    for pair, c in occ.most_common(5):
        print("  %s-%s  %.1f%%" % (pair[0], pair[1], 100.0 * c / n))
    np.savetxt(out_csv.replace(".csv", "_per_frame.txt"),
               np.column_stack([np.arange(n) * 0.3, counts]), fmt="%.1f")


analyse(r"D:\user data\Desktop\斑马鱼\5.dockMD\MD\03_production\traj167.pdb",
        "BDPrB", r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\hbond_BDPrB.csv")
analyse(r"D:\user data\Desktop\斑马鱼\5.dockMD\controls\kya33\traj167.pdb",
        "KyA33", r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\hbond_KyA33.csv")
