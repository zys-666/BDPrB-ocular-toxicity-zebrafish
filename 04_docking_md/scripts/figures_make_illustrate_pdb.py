# -*- coding: utf-8 -*-
"""Build combined PDBs (superposed / side-by-side) of the zebrafish model and human 9R6N
for the Illustrate web server (PDB only, no CIF)."""
import numpy as np, os, math

root = r"D:\user data\Desktop\斑马鱼\5.dockMD"
zpath = os.path.join(root, "dock", "zebrafish_ALDH1A2_trim26-518.pdb")
hpath = os.path.join(root, "dock", "recontrol", "protein.pdb")
out = os.path.join(root, "figures")
os.makedirs(out, exist_ok=True)

def read_pdb(path):
    atoms = []
    for l in open(path):
        if l.startswith("ATOM"):
            atoms.append(l.rstrip("\n"))
    return atoms

def ca(atoms, chain=None):
    d = {}
    for l in atoms:
        if l[12:16].strip() == "CA":
            d[int(l[22:26])] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return d

za = read_pdb(zpath)
ha = read_pdb(hpath)
zc = ca(za); hc = ca(ha)
OFF = 23  # 9R6N PDB numbering = zebrafish/UniProt numbering - 23
ks = sorted(k for k in hc if (k + OFF) in zc)
P = np.array([hc[k] for k in ks]); Q = np.array([zc[k + OFF] for k in ks])
Pm, Qm = P.mean(0), Q.mean(0)
Pc, Qc = P - Pm, Q - Qm
U, S, Vt = np.linalg.svd(Pc.T @ Qc)
d = np.sign(np.linalg.det(U @ Vt))
R = U @ np.diag([1, 1, d]) @ Vt
rms = math.sqrt((((Pc @ R) - Qc) ** 2).sum() / len(P))
print("CA superposition: %d residues, RMSD %.2f A" % (len(ks), rms))

def transform_line(l, R, Pm, Qm, shift=np.zeros(3)):
    x = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    y = (x - Pm) @ R + Qm + shift
    return l[:30] + "%8.3f%8.3f%8.3f" % (y[0], y[1], y[2]) + l[54:]

def write_pdb(path, groups):
    n = 0
    with open(path, "w") as f:
        for chain, lines, Rm, Pm, Qm, shift in groups:
            for l in lines:
                n += 1
                if Rm is None:
                    xyz = l[30:38] + l[38:46] + l[46:54]
                    body = l[:30] + xyz + l[54:]
                else:
                    body = transform_line(l, Rm, Pm, Qm, shift)
                f.write(body[:21] + chain + body[22:].replace(body[22:26], "%4d" % int(body[22:26]), 1))
                # rewrite serial properly
            f.write("TER\n")
        f.write("END\n")
    print("wrote", path, n, "atoms")

# simple serial renumbering version
def write_pdb2(path, blocks):
    n = 0
    with open(path, "w") as f:
        for chain, lines in blocks:
            for l in lines:
                n += 1
                l = l[:6] + "%5d" % n + l[11:]
                l = l[:21] + chain + l[22:]
                f.write(l + "\n")
            f.write("TER\n")
        f.write("END\n")
    print("wrote", path, n, "atoms")

z_lines = [l[:21] + "A" + l[22:] for l in za]
h_lines = [transform_line(l, R, Pm, Qm) for l in ha]
h_lines = [l[:21] + "B" + l[22:] for l in h_lines]

write_pdb2(os.path.join(out, "ill_zfish_vs_human_superposed.pdb"),
           [("A", z_lines), ("B", h_lines)])

shift = np.array([90.0, 0.0, 0.0])
h_side = [transform_line(l, R, Pm, Qm, shift) for l in ha]
h_side = [l[:21] + "B" + l[22:] for l in h_side]
write_pdb2(os.path.join(out, "ill_zfish_vs_human_sidebyside.pdb"),
           [("A", z_lines), ("B", h_side)])
