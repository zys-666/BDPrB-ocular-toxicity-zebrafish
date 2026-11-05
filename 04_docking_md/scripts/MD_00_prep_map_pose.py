# -*- coding: utf-8 -*-
"""Map the Vina pose back onto the MOL2 atom order using pre-dock coordinates as fingerprint.
Run from 00_prep/."""
import numpy as np, os

dock = os.path.join('..', '..', 'dock', 'dock_run')

def coords(path):
    out = []
    for l in open(path):
        if l.startswith(('ATOM', 'HETATM')):
            out.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    return np.array(out)

orig = coords(os.path.join(dock, 'ligand.pdbqt'))
print('pre-dock atoms:', orig.shape)

models, cur = [], None
for l in open(os.path.join(dock, 'dock_out.pdbqt')):
    if l.startswith('MODEL'):
        cur = []
    elif l.startswith('ENDMDL'):
        models.append(np.array(cur)); cur = None
    elif l.startswith(('ATOM', 'HETATM')) and cur is not None:
        cur.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
pose = models[0]
print('pose atoms:', pose.shape)

lines = open(os.path.join(dock, 'ligand.mol2')).read().splitlines()
i = lines.index('@<TRIPOS>ATOM'); j = lines.index('@<TRIPOS>BOND')
atoms = [l.split() for l in lines[i+1:j] if l.strip()]
bonds = [l.split() for l in lines[j+1:] if l.strip() and not l.startswith('@')]
mol2_heavy = np.array([[float(p[2]), float(p[3]), float(p[4])] for p in atoms[:22]])

perm = []
for k in range(22):
    d = np.linalg.norm(mol2_heavy - orig[k], axis=1)
    m = int(np.argmin(d))
    assert d[m] < 0.02, (k, m, d[m])
    perm.append(m)
print('permutation (pdbqt -> mol2):', perm)
assert sorted(perm) == list(range(22)), 'not a bijection'

newcoords = np.zeros((22, 3))
for k in range(22):
    newcoords[perm[k]] = pose[k]

heavy_bonds = [b for b in bonds if int(b[1]) <= 22 and int(b[2]) <= 22]
out = ['@<TRIPOS>MOLECULE', 'BDPrB_docked', '%4d %4d 0 0 0' % (22, len(heavy_bonds)),
       'SMALL', 'USER_CHARGES', '', '@<TRIPOS>ATOM']
for a in range(22):
    p = atoms[a]
    x, y, z = newcoords[a]
    out.append('%7d %-6s %9.4f %9.4f %9.4f %-6s 1 UNL1 %9.4f' % (a+1, p[1], x, y, z, p[5], float(p[8])))
out.append('@<TRIPOS>BOND')
for b in heavy_bonds:
    out.append('%6s %4d %4d %-4s' % (b[0], int(b[1]), int(b[2]), b[3]))
out.append('')
open('lig_mapped_heavy.mol2', 'w').write('\n'.join(out) + '\n')
print('wrote lig_mapped_heavy.mol2 with', len(heavy_bonds), 'bonds')

for b in heavy_bonds[:6]:
    a1, a2 = int(b[1])-1, int(b[2])-1
    print('  bond %2d-%2d %-4s %-4s len %.2f' % (a1+1, a2+1, atoms[a1][5], atoms[a2][5],
                                                 np.linalg.norm(newcoords[a1]-newcoords[a2])))
