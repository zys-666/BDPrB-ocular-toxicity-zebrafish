# -*- coding: utf-8 -*-
"""Symmetry-corrected heavy-atom RMSD between the redocked KyA33 pose and the crystal pose."""
import numpy as np, math, sys

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\dock\recontrol"

def read_mol2(path):
    lines = open(path).read().splitlines()
    i = lines.index('@<TRIPOS>ATOM'); j = lines.index('@<TRIPOS>BOND')
    atoms = []
    for l in lines[i+1:j]:
        p = l.split()
        if len(p) >= 6:
            atoms.append((p[5].split('.')[0], np.array([float(p[2]), float(p[3]), float(p[4])])))
    bonds = []
    for l in lines[j+1:]:
        p = l.split()
        if len(p) >= 4 and not l.startswith('@'):
            bonds.append((int(p[1]) - 1, int(p[2]) - 1))
    return atoms, bonds

def build(path):
    at, bo = read_mol2(path)
    heavy = [k for k, a in enumerate(at) if a[0] != 'H']
    idx = {k: n for n, k in enumerate(heavy)}
    X = np.array([at[k][1] for k in heavy])
    el = [at[k][0] for k in heavy]
    adj = [set() for _ in heavy]
    for a, b in bo:
        if a in idx and b in idx:
            adj[idx[a]].add(idx[b]); adj[idx[b]].add(idx[a])
    return X, el, adj

Xc, elc, adjc = build(base + r"\lig_ref.mol2")
Xd, eld, adjd = build(base + r"\pose1_noh.mol2")
print("crystal: %d heavy atoms, docked: %d" % (len(Xc), len(Xd)))
assert len(Xc) == len(Xd)

deg = [len(a) for a in adjc]
degd = [len(a) for a in adjd]

# backtracking graph isomorphism (connectivity + element + degree)
n = len(Xc)
order = sorted(range(n), key=lambda i: (-deg[i], elc[i]))
best = [None]
best_rmsd = [1e9]

def kabsch_rmsd(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(U @ Vt))
    R = U @ np.diag([1, 1, d]) @ Vt
    return math.sqrt((((Pc @ R) - Qc) ** 2).sum() / len(P))

def search(k, mapping):
    if best_rmsd[0] < 0.05:
        return
    if k == n:
        P = np.array([Xc[i] for i in range(n)])
        Q = np.array([Xd[mapping[i]] for i in range(n)])
        r = kabsch_rmsd(P, Q)
        if r < best_rmsd[0]:
            best_rmsd[0] = r
            best[0] = dict(mapping)
        return
    i = order[k]
    used = set(mapping.values())
    for j in range(n):
        if j in used: continue
        if eld[j] != elc[i] or degd[j] != deg[i]: continue
        ok = True
        for nb in adjc[i]:
            if nb in mapping:
                if mapping[nb] not in adjd[j]:
                    ok = False; break
        if not ok: continue
        mapping[i] = j
        search(k + 1, mapping)
        del mapping[i]

search(0, {})
print("best symmetry-corrected heavy-atom RMSD: %.2f A" % best_rmsd[0])
mp = best[0]
# centroid distance
P = np.array([Xc[i] for i in range(n)])
Q = np.array([Xd[mp[i]] for i in range(n)])
print("crystal centroid:", P.mean(0).round(2), " docked centroid:", Q.mean(0).round(2))
print("centroid distance: %.2f A" % np.linalg.norm(P.mean(0) - Q.mean(0)))
