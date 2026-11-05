# -*- coding: utf-8 -*-
"""Merge protein (pdb2gmx) and ligand (acpype) topologies + coordinates."""
import os, re

base = os.path.dirname(os.path.abspath(__file__))
itp_path = os.path.join(base, 'LIG.itp')
lines = open(itp_path).read().splitlines()

# split atomtypes block
at_start = next(i for i, l in enumerate(lines) if l.strip() == '[ atomtypes ]')
mol_start = next(i for i, l in enumerate(lines) if l.strip() == '[ moleculetype ]')

atomtypes = lines[at_start:mol_start]
molecule = lines[mol_start:]

# rename residue UNL -> LIG inside the molecule block (atom lines only)
newmol = []
for l in molecule:
    if re.match(r'^\s+\d+\s+\S+\s+1\s+UNL', l):
        l = l.replace('UNL', 'LIG')
    newmol.append(l)

open(os.path.join(base, 'LIG_atomtypes.itp'), 'w').write('\n'.join(atomtypes) + '\n')
open(os.path.join(base, 'LIG.itp'), 'w').write('\n'.join(newmol) + '\n')
print('split: atomtypes lines %d, molecule lines %d' % (len(atomtypes), len(newmol)))

# ---- complex.gro ----
def read_gro(path):
    L = open(path).read().splitlines()
    n = int(L[1])
    atoms = L[2:2+n]
    box = L[2+n]
    return n, atoms, box

np_, ap, bp = read_gro(os.path.join(base, 'protein.gro'))
nl, al, bl = read_gro(os.path.join(base, 'ligand', 'LIG_GMX.gro'))
# rewrite ligand residue name UNL -> LIG
al2 = []
for l in al:
    al2.append('%5d%-5s%5s%5d' % (1, 'LIG', l[10:15].strip(), int(l[15:20])) + l[20:])
out = ['BDPrB-ALDH1A2 complex', '%5d' % (np_ + nl)] + ap + al2 + [bp]
open(os.path.join(base, 'complex.gro'), 'w').write('\n'.join(out) + '\n')
print('complex.gro: %d + %d = %d atoms, box %s' % (np_, nl, np_ + nl, bp))
