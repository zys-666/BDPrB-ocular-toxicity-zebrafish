# -*- coding: utf-8 -*-
"""Build the protein-ligand complex input for MD from the Vina pose."""
import os, numpy as np

dock = r"D:\user data\Desktop\斑马鱼\5.dockMD\dock"
md = r"D:\user data\Desktop\斑马鱼\5.dockMD\MD"
prep = os.path.join(md, '00_prep')

# ---- 1. ligand pose (model 1) from Vina ----
models = []
cur = None
for l in open(os.path.join(dock, 'dock_run', 'dock_out.pdbqt')):
    if l.startswith('MODEL'):
        cur = []
    elif l.startswith('ENDMDL'):
        models.append(cur); cur = None
    elif l.startswith(('ATOM', 'HETATM')) and cur is not None:
        cur.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
pose = np.array(models[0])
print('pose atoms:', pose.shape)

# ---- 2. mol2 template (topology) ----
lines = open(os.path.join(dock, 'dock_run', 'ligand.mol2')).read().splitlines()
mol_start = lines.index('@<TRIPOS>MOLECULE')
atom_start = lines.index('@<TRIPOS>ATOM')
bond_start = lines.index('@<TRIPOS>BOND')
atoms = [l.split() for l in lines[atom_start+1:bond_start] if l.strip()]
bonds = [l.split() for l in lines[bond_start+1:] if l.strip() and not l.startswith('@')]
print('mol2 atoms:', len(atoms), 'bonds:', len(bonds))
assert len(atoms) == 22 + 28

out = []
out.append('@<TRIPOS>MOLECULE')
out.append('BDPrB_docked_pose')
nb_written = sum(1 for b in bonds if int(b[1])<=22 and int(b[2])<=22)
out.append('%4d %4d 0 0 0' % (22, nb_written))
print('heavy bonds:', nb_written)
out.append('SMALL')
out.append('USER_CHARGES')
out.append('')
out.append('@<TRIPOS>ATOM')
for i in range(22):
    p = atoms[i]
    x, y, z = pose[i]
    out.append('%7d %-6s %9.4f %9.4f %9.4f %-6s 1 UNL1 %9.4f' % (
        i+1, p[1], x, y, z, p[5], float(p[8])))
out.append('@<TRIPOS>BOND')
for b in bonds:
    i, j = int(b[1]), int(b[2])
    if i <= 22 and j <= 22:
        out.append('%6s %4d %4d %-4s' % (b[0], i, j, b[3]))
out.append('')
open(os.path.join(prep, 'lig_pose_heavy.mol2'), 'w').write('\n'.join(out) + '\n')
print('wrote lig_pose_heavy.mol2')

# ---- 3. complex pdb (protein + ligand heavy atoms) ----
prot = [l for l in open(os.path.join(dock, 'dock_run', 'protein.pdb')) if l.startswith(('ATOM', 'HETATM'))]
with open(os.path.join(prep, 'complex.pdb'), 'w') as f:
    for l in prot:
        f.write(l)
    f.write('TER\n')
    for i in range(22):
        x, y, z = pose[i]
        el = atoms[i][5].split('.')[0][0]
        f.write('HETATM%5d %-4s LIG L 900    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n' % (
            i+1, atoms[i][1], x, y, z, el))
    f.write('END\n')
print('wrote complex.pdb with', len(prot), 'protein atoms + 22 ligand atoms')

# centroid check
lig_c = pose.mean(0)
print('ligand centroid in complex:', lig_c.round(2))
