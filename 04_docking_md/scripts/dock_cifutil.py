# -*- coding: utf-8 -*-
"""Shared mmCIF parsing / alignment / superposition helpers."""
# -*- coding: utf-8 -*-
"""Compare human ALDH1A2/KyA33 crystal structure (9R6N) with zebrafish ALDH1A2 models."""
import sys, math, re, json
import numpy as np

AA3 = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}

def tokenize(line):
    out=[]; i=0; n=len(line)
    while i<n:
        c=line[i]
        if c in ' \t':
            i+=1; continue
        if c in "'\"":
            q=c; i+=1; j=line.find(q,i)
            out.append(line[i:j]); i=j+1; continue
        j=i
        while j<n and line[j] not in ' \t':
            j+=1
        out.append(line[i:j]); i=j
    return out

def parse_cif(path):
    """Return dict of category -> list of dict rows (for loop_ cats) and scalar items."""
    txt = open(path, encoding='utf-8', errors='replace').read()
    lines = txt.splitlines()
    cats = {}
    scalars = {}
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        s = line.strip()
        if s.startswith('loop_'):
            i += 1
            heads = []
            while i < n and lines[i].strip().startswith('_'):
                heads.append(lines[i].strip())
                i += 1
            rows = []
            while i < n:
                s2 = lines[i].strip()
                if s2.startswith('#') or s2.startswith('loop_') or s2.startswith('_') or s2.startswith('data_'):
                    break
                if s2 == '' or s2.startswith(';'):
                    i += 1
                    continue
                vals = tokenize(lines[i])
                # handle semicolon text blocks crudely
                rows.append(vals)
                i += 1
            cat = heads[0].split('.')[0]
            cats.setdefault(cat, []).extend(rows)
            cats[cat + '__heads'] = heads
            continue
        if s.startswith('_'):
            parts = tokenize(s)
            key = parts[0]
            if len(parts) > 1:
                scalars[key] = ' '.join(parts[1:])
            else:
                # value on next line
                if i+1 < n:
                    scalars[key] = lines[i+1].strip()
                    i += 1
            i += 1
            continue
        i += 1
    return cats, scalars

def atoms(path):
    cats, scal = parse_cif(path)
    heads = cats['_atom_site__heads']
    idx = {h.split('.')[1]: k for k, h in enumerate(heads)}
    out = []
    for r in cats['_atom_site']:
        if len(r) < len(heads):
            r = r + ['?']*(len(heads)-len(r))
        try:
            out.append(dict(
                group=r[idx['group_PDB']],
                atom=r[idx['label_atom_id']].strip('"'),
                comp=r[idx['label_comp_id']],
                asym=r[idx['label_asym_id']],
                seq=r[idx['label_seq_id']],
                auth_seq=r[idx['auth_seq_id']],
                x=float(r[idx['Cartn_x']]), y=float(r[idx['Cartn_y']]), z=float(r[idx['Cartn_z']]),
                b=float(r[idx['B_iso_or_equiv']]) if r[idx['B_iso_or_equiv']] not in ('?','.') else float('nan'),
                alt=r[idx['label_alt_id']] if 'label_alt_id' in idx else '.',
            ))
        except Exception as e:
            continue
    return out, cats, scal

def polymer_residues(at):
    """ordered unique (asym, seq) -> dict with resname, ca coord, b"""
    res = {}
    order = []
    for a in at:
        if a['group'] != 'ATOM':
            continue
        key = (a['asym'], a['seq'])
        if key not in res:
            aa = AA3.get(a['comp'])
            if aa is None:
                continue
            res[key] = {'name': a['comp'], 'aa': aa, 'ca': None, 'b': a['b'], 'atoms': []}
            order.append(key)
        res[key]['atoms'].append(a)
        if a['atom'] == 'CA':
            res[key]['ca'] = np.array([a['x'], a['y'], a['z']])
            res[key]['b'] = a['b']
    return [(k, res[k]) for k in order if k in res]

def nw_align(s1, s2, match=2, mismatch=-1, gap=-7):
    n, m = len(s1), len(s2)
    NEG = -10**9
    # band-less full DP, 518x518 fine
    H = np.full((n+1, m+1), NEG, dtype=np.int32)
    P = np.zeros((n+1, m+1), dtype=np.int8)  # 0 diag 1 up 2 left
    H[0, :] = np.arange(m+1)*gap
    H[:, 0] = np.arange(n+1)*gap
    P[0, 1:] = 2
    P[1:, 0] = 1
    for i in range(1, n+1):
        c1 = s1[i-1]
        row = H[i]; prev = H[i-1]; prow = P[i]
        for j in range(1, m+1):
            d = prev[j-1] + (match if c1 == s2[j-1] else mismatch)
            u = prev[j] + gap
            l = row[j-1] + gap
            if d >= u and d >= l:
                row[j] = d; prow[j] = 0
            elif u >= l:
                row[j] = u; prow[j] = 1
            else:
                row[j] = l; prow[j] = 2
    a1, a2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        p = P[i, j]
        if i > 0 and j > 0 and p == 0:
            a1.append(s1[i-1]); a2.append(s2[j-1]); i -= 1; j -= 1
        elif i > 0 and p == 1:
            a1.append(s1[i-1]); a2.append('-'); i -= 1
        else:
            a1.append('-'); a2.append(s2[j-1]); j -= 1
    return ''.join(reversed(a1)), ''.join(reversed(a2))

def kabsch(P, Q):
    Pc = P - P.mean(0); Qc = Q - Q.mean(0)
    C = Pc.T @ Qc
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1, 1, d])
    R = V @ D @ Wt
    rmsd = math.sqrt((((Pc @ R) - Qc)**2).sum() / len(P))
    return R, P.mean(0), Q.mean(0), rmsd


