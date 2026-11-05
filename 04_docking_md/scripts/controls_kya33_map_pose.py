# -*- coding: utf-8 -*-
"""Map the docked KyA33 pose (zebrafish) back onto the mol2 atom order."""
import numpy as np, os, sys

dock = r"D:\user data\Desktop\斑马鱼\5.dockMD\dock"
out = r"D:\user data\Desktop\斑马鱼\5.dockMD\controls\kya33"
os.makedirs(out, exist_ok=True)

def coords(path):
    r = []
    for l in open(path):
        if l.startswith(('ATOM', 'HETATM')):
            el = l[77:79].strip()
            if el in ('H', 'HD'):
                continue
            r.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(r)

orig = coords(os.path.join(dock, "recontrol", "kya33.pdbqt"))
pose = coords(os.path.join(dock, "dock_run", "kya33_zebrafish.pdbqt"))
# first model only
n_at = len(orig)
pose = pose[:n_at]
print("atoms:", orig.shape, pose.shape)

lines = open(os.path.join(dock, "recontrol", "kya33.mol2")).read().splitlines()
i = lines.index('@<TRIPOS>ATOM'); j = lines.index('@<TRIPOS>BOND')
atoms = [l.split() for l in lines[i+1:j] if l.strip()]
bonds = [l.split() for l in lines[j+1:] if l.strip() and not l.startswith('@')]
heavy = [k for k, a in enumerate(atoms) if a[5].split('.')[0] != 'H']
print("mol2 atoms:", len(atoms), "heavy:", len(heavy))

mol2_heavy = np.array([[float(atoms[k][2]), float(atoms[k][3]), float(atoms[k][4])] for k in heavy])
perm = []
for a in range(n_at):
    d = np.linalg.norm(mol2_heavy - orig[a], axis=1)
    m = int(np.argmin(d))
    assert d[m] < 0.05, (a, m, d[m])
    perm.append(heavy[m])
print("bijection:", sorted(perm) == sorted(heavy))

newc = np.zeros((len(atoms), 3))
for a in range(n_at):
    newc[perm[a]] = pose[a]

heavy_bonds = [b for b in bonds if int(b[1]) - 1 in heavy and int(b[2]) - 1 in heavy]
outl = ['@<TRIPOS>MOLECULE', 'KYA33_docked', '%4d %4d 0 0 0' % (len(heavy), len(heavy_bonds)),
        'SMALL', 'USER_CHARGES', '', '@<TRIPOS>ATOM']
for n, k in enumerate(heavy, 1):
    a = atoms[k]
    outl.append('%7d %-6s %9.4f %9.4f %9.4f %-6s 1 UNL1 %9.4f' % (
        n, a[1], newc[k][0], newc[k][1], newc[k][2], a[5], float(a[8])))
outl.append('@<TRIPOS>BOND')
idmap = {k: n for n, k in enumerate(heavy, 1)}
for b in heavy_bonds:
    outl.append('%6s %4d %4d %-4s' % (b[0], idmap[int(b[1]) - 1], idmap[int(b[2]) - 1], b[3]))
outl.append('')
open(os.path.join(out, 'kya33_heavy.mol2'), 'w').write('\n'.join(outl) + '\n')
print("wrote kya33_heavy.mol2 with", len(heavy), "atoms")
