# -*- coding: utf-8 -*-
"""Build the MM-PBSA index for the KyA33 control (Protein / LIGAND / Protein_LIGAND)."""
import os

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\controls\kya33"

def read_ndx(path):
    groups = {}
    name = None
    for line in open(path):
        s = line.strip()
        if s.startswith('[') and s.endswith(']'):
            name = s[1:-1].strip()
            groups[name] = []
        elif name and s:
            groups[name].extend(int(x) for x in s.split())
    return groups

defg = read_ndx(os.path.join(base, 'def.ndx'))
lig = read_ndx(os.path.join(base, 'lig.ndx'))
prot = defg['Protein']
ligand = list(lig.values())[0]
combo = sorted(set(prot) | set(ligand))

def write(f, nm, atoms):
    f.write('[ %s ]\n' % nm)
    for i in range(0, len(atoms), 15):
        f.write(' '.join(str(a) for a in atoms[i:i+15]) + '\n')

with open(os.path.join(base, 'index.ndx'), 'w') as f:
    write(f, 'Protein', prot)
    write(f, 'LIGAND', ligand)
    write(f, 'Protein_LIGAND', combo)
    write(f, 'System', defg['System'])
print('protein %d, ligand %d, complex %d' % (len(prot), len(ligand), len(combo)))
