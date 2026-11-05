# -*- coding: utf-8 -*-
"""Prepare a docking-ready zebrafish ALDH1A2 receptor from the AF model + reference ligand pose."""
import os, math
import numpy as np
from cifutil import atoms, polymer_residues, nw_align, kabsch

base = os.path.dirname(os.path.abspath(__file__))
h_at, _, _ = atoms(os.path.join(base, '9R6N.cif'))
h_res = polymer_residues(h_at)
z_at, _, _ = atoms(os.path.join(base, 'AF-Q90XS8-F1-model_v6.cif'))
z_res = polymer_residues(z_at)
hA = {int(s): r for (a, s), r in h_res if a == 'A'}
zA = {int(s): r for (a, s), r in z_res if a == 'A'}
hk = sorted(hA); zk = sorted(zA)
a1, a2 = nw_align(''.join(hA[k]['aa'] for k in hk), ''.join(zA[k]['aa'] for k in zk))
h2z = {}; ih = iz = 0
for c1, c2 in zip(a1, a2):
    if c1 != '-': ih += 1
    if c2 != '-': iz += 1
    if c1 != '-' and c2 != '-': h2z[hk[ih-1]] = zk[iz-1]

P = np.array([zA[h2z[k]]['ca'] for k in hk]); Q = np.array([hA[k]['ca'] for k in hk])
R, cP, cQ, rmsd = kabsch(P, Q)
t = cQ - cP @ R
lig = np.array([[a['x'], a['y'], a['z']] for a in h_at
                if a['group'] == 'HETATM' and a['comp'] == 'A1IDQ' and a['asym'] == 'O'])
lig_z = (lig - t) @ R.T
cen = lig_z.mean(0)

TRIM = 25
keep = [k for k in zk if k > TRIM]
Z = []
lbl = []
for k in keep:
    for a in zA[k]['atoms']:
        if a['atom'] in ('N', 'CA', 'C', 'O'): continue
        Z.append([a['x'], a['y'], a['z']]); lbl.append((k, a['atom'], a['comp']))
Z = np.array(Z)
D = np.sqrt(((lig_z[:, None, :] - Z[None, :, :])**2).sum(-1))
print('after trimming 1-%d: min dist ligand-heavy to side chain = %.2f A' % (TRIM, D.min()))
print('clashes < 2.2 A:', int((D.min(0) < 2.2).sum()))
print('contacts < 3.6 A:')
seen = set()
for i in np.argsort(D.min(0)):
    d = D[:, i].min()
    if d > 3.6: break
    k, an, cn = lbl[i]
    if (k, cn) in seen: continue
    seen.add((k, cn))
    print('   Z%3d %s %-4s %.2f' % (k, cn, an, d))

print('\nsubstituted pocket positions (zebrafish numbering):')
for k in (191, 477, 478, 479):
    ds = [D[:, i].min() for i, (k2, an, cn) in enumerate(lbl) if k2 == k]
    print('   Z%3d %s  min dist to KyA33 = %.2f  pLDDT %.1f' % (k, zA[k]['name'], min(ds), zA[k]['b']))

# write trimmed receptor PDB (no hydrogens, no ligand)
out = []
atomserial = 1
for k in keep:
    for a in zA[k]['atoms']:
        el = a['atom'][0] if a['atom'][0].isalpha() else 'C'
        out.append('ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f           %2s' % (
            atomserial, a['atom'], a['comp'], k, a['x'], a['y'], a['z'], a['b'], el))
        atomserial += 1
out.append('END')
open(os.path.join(base, 'zebrafish_ALDH1A2_trim26-518.pdb'), 'w').write('\n'.join(out) + '\n')
print('\nwrote zebrafish_ALDH1A2_trim26-518.pdb with %d atoms, %d residues' % (atomserial-1, len(keep)))

# write transplanted reference ligand pose
lines = []
for i, (x, y, z) in enumerate(lig_z, 1):
    lines.append('HETATM%5d  C   LIG A 900    %8.3f%8.3f%8.3f  1.00  0.00           C' % (i, x, y, z))
lines.append('END')
open(os.path.join(base, 'KyA33_reference_pose_zebrafish.pdb'), 'w').write('\n'.join(lines) + '\n')
print('wrote KyA33_reference_pose_zebrafish.pdb')
print('pocket centre (zebrafish frame): %.2f %.2f %.2f' % tuple(cen))
print('ligand centroid radial extent: %.2f A' % np.linalg.norm(lig_z - cen, axis=1).max())

# sphere of pocket residues: max CA distance from centre
ds = [(k, np.linalg.norm(zA[k]['ca'] - cen)) for k in keep]
ds.sort(key=lambda x: x[1])
print('nearest residue CAs:', [(k, round(d, 1)) for k, d in ds[:6]])
