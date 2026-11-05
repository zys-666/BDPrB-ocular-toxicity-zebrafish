# -*- coding: utf-8 -*-
"""Extract chain A protein and its KyA33 ligand from 9R6N for a redocking control."""
import numpy as np, os, sys

sys.path.insert(0, r"D:\user data\Desktop\斑马鱼\5.dockMD\dock")
from cifutil import atoms

src = r"D:\user data\Desktop\斑马鱼\5.dockMD\dock\9R6N.cif"
out = r"D:\user data\Desktop\斑马鱼\5.dockMD\dock\recontrol"
os.makedirs(out, exist_ok=True)

at = atoms(src)[0]
prot = [a for a in at if a['group'] == 'ATOM' and a['asym'] == 'A']
lig = [a for a in at if a['group'] == 'HETATM' and a['comp'] == 'A1IDQ']
# ligand copy nearest chain A
ligpts = np.array([[a['x'], a['y'], a['z']] for a in lig])
protpts = np.array([[a['x'], a['y'], a['z']] for a in prot])
best = None
for copy in sorted(set(a['asym'] for a in lig)):
    L = np.array([[a['x'], a['y'], a['z']] for a in lig if a['asym'] == copy])
    d = np.linalg.norm(L[:, None] - protpts[None], axis=-1).min()
    if best is None or d < best[1]:
        best = (copy, d)
print('ligand copy nearest chain A:', best)
L = [a for a in lig if a['asym'] == best[0]]

with open(os.path.join(out, 'protein.pdb'), 'w') as f:
    for i, a in enumerate(prot, 1):
        f.write('ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s\n' % (
            i, a['atom'], a['comp'], int(a['seq']), a['x'], a['y'], a['z'], a['b'],
            (a['atom'][0] if a['atom'][0].isalpha() else 'C')))
    f.write('END\n')
with open(os.path.join(out, 'lig_ref.pdb'), 'w') as f:
    for i, a in enumerate(L, 1):
        f.write('HETATM%5d %-4s LIG L 900    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n' % (
            i, a['atom'], a['x'], a['y'], a['z'],
            (a['atom'][0] if a['atom'][0].isalpha() else 'C')))
    f.write('END\n')
cen = np.array([[a['x'], a['y'], a['z']] for a in L]).mean(0)
print('protein atoms:', len(prot), ' ligand atoms:', len(L), ' ligand centroid:', cen.round(2))
print('ligand extent:', (np.array([[a['x'], a['y'], a['z']] for a in L]).max(0) -
                         np.array([[a['x'], a['y'], a['z']] for a in L]).min(0)).round(1))
