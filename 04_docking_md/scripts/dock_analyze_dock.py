# -*- coding: utf-8 -*-
"""Analyse Vina output poses: centroid, pocket contacts, distance to Cys320 SG."""
import numpy as np, math, sys, os, glob

base = os.path.dirname(os.path.abspath(__file__))

def read_pdbqt_models(path):
    models = []
    cur = None
    for l in open(path):
        if l.startswith('MODEL'):
            cur = {'atoms': [], 'score': None, 'idx': int(l.split()[1])}
        elif l.startswith('ENDMDL'):
            models.append(cur); cur = None
        elif l.startswith(('ATOM', 'HETATM')) and cur is not None:
            cur['atoms'].append((l[12:16].strip(), l[17:20].strip(), int(l[22:26]),
                                 np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]),
                                 l[77:79].strip()))
        elif l.startswith('REMARK VINA RESULT') and cur is not None:
            cur['score'] = float(l.split()[3])
    return models

def pocket_center():
    L = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                  for l in open(os.path.join(base, 'KyA33_reference_pose_zebrafish.pdb')) if l.startswith('HETATM')])
    return L.mean(0), L

# receptor atoms of interest
rec = []
for l in open(os.path.join(base, 'dock_run', 'protein.pdb')):
    if l.startswith('ATOM'):
        rec.append((int(l[22:26]), l[12:16].strip(), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))

def min_dist_to(resnum, atomname, X):
    sel = [x for (r, n, x) in rec if r == resnum and n == atomname]
    if not sel:
        return float('nan')
    p = sel[0]
    return np.min(np.linalg.norm(X - p, axis=1))

def contacts(X, cut=4.0):
    out = {}
    for (r, n, x) in rec:
        d = np.min(np.linalg.norm(X - x, axis=1))
        if d <= cut:
            out.setdefault(r, []).append((n, d))
    return out

center, Lref = pocket_center()
for f in sorted(glob.glob(os.path.join(base, 'dock_run', 'dock_out*.pdbqt'))):
    models = read_pdbqt_models(f)
    print('\n==== %s : %d models ====' % (os.path.basename(f), len(models)))
    print('%4s %8s %8s %9s %9s %6s %6s' % ('mode', 'score', 'c2pocket', 'd(C320SG)', 'd(C319SG)', 'ncont', 'nH'))
    for m in models[:10]:
        X = np.array([a[3] for a in m['atoms']])
        heavy = X[[i for i, a in enumerate(m['atoms']) if a[4] != 'HD']]
        c = heavy.mean(0)
        cont = contacts(heavy, 4.0)
        print('%4d %8.1f %8.2f %9.2f %9.2f %6d' % (
            m['idx'], m['score'], np.linalg.norm(c - center),
            min_dist_to(320, 'SG', heavy), min_dist_to(319, 'SG', heavy), len(cont)))
    # detailed contacts for best
    m = models[0]
    heavy = np.array([a[3] for a in m['atoms'] if a[4] != 'HD'])
    cont = contacts(heavy, 3.8)
    print('  best pose contacts <3.8 A:')
    for r in sorted(cont):
        atoms = sorted(set(n for n, d in cont[r]))
        dmin = min(d for n, d in cont[r])
        nm = [x[1] for x in rec if x[0] == r][0]
        print('    %3d %-4s %-30s %.2f' % (r, nm, ','.join(atoms[:6]), dmin))
