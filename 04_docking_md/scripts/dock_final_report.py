# -*- coding: utf-8 -*-
"""Final comparison: 9R6N (human ALDH1A2 + KyA33) vs zebrafish ALDH1A2 AF models."""
import os, sys, math, json
import numpy as np
from cifutil import atoms, polymer_residues, nw_align, kabsch

base = os.path.dirname(os.path.abspath(__file__))
human_uniprot = ("MTSSKIEMPGEVKADPAALMASLHLLPSPTPNLEIKYTKIFINNEWQNSESGRVFPVYNPATGEQVCEVQEADKADIDKAVQAARLAFSLGSVWRRMDASERGRLLDKLADLVERDRAVL"
"ATMESLNGGKPFLQAFYVDLQGVIKTFRYYAGWADKIHGMTIPVDGDYFTFTRHEPIGVCGQIIPWNFPLLMFAWKIAPALCCGNTVVIKPAEQTPLSALYMGALIKEAGFPPGVINILP"
"GYGPTAGAAIASHIGIDKIAFTGSTEVGKLIQEAAGRSNLKRVTLELGGKSPNIIFADADLDYAVEQAHQGVFFNQGQCCTAGSRIFVEESIYEEFVRRSVERAKRRVVGSPFDPTTEQG"
"PQIDKKQYNKILELIQSGVAEGAKLECGGKGLGRKGFFIEPTVFSNVTDDMRIAKEEIFGPVQEILRFKTMDEVIERANNSDFGLVAAVFTNDINKALTVSSAMQAGTVWINCYNALNAQ"
"SPFGGFKMSGNGREMGEFGLREYSEVKTVTVKIPQKNS")
zfish_uniprot = ("MTSSEVELPGEVKNDPAALMASLHLMPSPVPNPEIKYTKIFINNEWHDSVSGKVFPTYNPATGEKICDVQEADKADVDKAVQAARSAFSLGSVWRKMDASERGKLLFKLADLVERDSAYL"
"ATLESLDSGKPFLPCFFVDLQGIIKTFRYYAGWADKIHGSTIPIDGEFFTLTRHEPIGVCGQIIPWNFPLVMTAWKLGPALSCGNTVVLKPAEQTPLTCLYLGALIKEAGFPPGVVNILP"
"GYGPTAGAAISSHMGIDKVAFTGSTEVGKLIQEAAGKSNLKRVTLELGGKSPNIIFADADFELALEQAHQGVFFNNGQCCTAGSRIFVEEPIYDEFVRRSVERAQRRKVGNPFDPTTEHG"
"PQVSEEQQRRVLELIQSGITEGAKLECGGKAPATKGFFVEPTVFSNVKDHMRIAKEEIFGPVQQIMKFKTIEEVIERANNTEYGLAAAVFTRDISKAMTISAAVQAGTVWINCYNAMSCQ"
"CPFGGFKMSGNGRELGEIGLKEYTELKTITMKMSGKTS")
assert len(human_uniprot) == 518, len(human_uniprot)
assert len(zfish_uniprot) == 518, len(zfish_uniprot)

h_at, _, _ = atoms(os.path.join(base, '9R6N.cif'))
h_res = polymer_residues(h_at)
z_at, _, _ = atoms(os.path.join(base, 'AF-Q90XS8-F1-model_v6.cif'))
z_res = polymer_residues(z_at)
z3_at, _, _ = atoms(os.path.join(base, 'AF-0000000066349646-model_v1.cif'))
z3_res = polymer_residues(z3_at)

hA = {int(s): r for (a, s), r in h_res if a == 'A'}
zA = {int(s): r for (a, s), r in z_res if a == 'A'}
z3A = {int(s): r for (a, s), r in z3_res if a == 'A'}
hk = sorted(hA); zk = sorted(zA)

# 1) map 9R6N numbering -> UniProt human numbering
a1, a2 = nw_align(''.join(hA[k]['aa'] for k in hk), human_uniprot)
offsets = {}
ih = 0; iu = 0
for c1, c2 in zip(a1, a2):
    if c1 != '-': ih += 1
    if c2 != '-': iu += 1
    if c1 != '-' and c2 != '-':
        offsets[hk[ih-1]] = iu
print('9R6N numbering -> UniProt offset sample:', [(k, offsets[k]) for k in hk[:3]], [(k, offsets[k]) for k in hk[-3:]])
off = set(offsets[k] - k for k in hk)
print('constant offset?', off)

# 2) 9R6N -> zebrafish numbering
a1, a2 = nw_align(''.join(hA[k]['aa'] for k in hk), ''.join(zA[k]['aa'] for k in zk))
h2z = {}
ih = iz = 0
for c1, c2 in zip(a1, a2):
    if c1 != '-': ih += 1
    if c2 != '-': iz += 1
    if c1 != '-' and c2 != '-':
        h2z[hk[ih-1]] = zk[iz-1]
print('pairs', len(h2z))

nid = sum(1 for x, y in zip(a1, a2) if x == y and x != '-')
print('zebrafish vs 9R6N identity: %d/%d = %.1f%%' % (nid, len(h2z), 100*nid/len(h2z)))
# full-length uniprot identity
b1, b2 = nw_align(human_uniprot, zfish_uniprot)
ni = sum(1 for x, y in zip(b1, b2) if x == y and x != '-')
print('full-length human vs zebrafish(UniProt): %d/%d = %.1f%%, gaps %d/%d' % (ni, len(human_uniprot), 100*ni/len(human_uniprot), b1.count('-'), b2.count('-')))

# 3) pocket definition
lig = np.array([[a['x'], a['y'], a['z']] for a in h_at if a['group'] == 'HETATM' and a['comp'] == 'A1IDQ' and a['asym'] == 'O'])
# ligand chain O nearest to chain A
protA = np.array([[q['x'], q['y'], q['z']] for k in hk for q in hA[k]['atoms']])
d = np.sqrt(((protA[:, None, :] - lig[None])**2).sum(-1)).min()
print('ligand copy O min dist to chain A: %.2f' % d)

pk = []
for k in hk:
    ats = np.array([[q['x'], q['y'], q['z']] for q in hA[k]['atoms']])
    dd = np.sqrt(((ats[:, None, :] - lig[None])**2).sum(-1)).min()
    pk.append((k, dd))
pocket = [(k, dd) for k, dd in pk if dd <= 4.5]
print('pocket size', len(pocket))

# 4) second shell 5-8 A
shell2 = [(k, dd) for k, dd in pk if 4.5 < dd <= 8.0]
id2 = sum(1 for k, dd in shell2 if zA[h2z[k]]['aa'] == hA[k]['aa'])
print('second shell (4.5-8 A): %d residues, identity %d/%d = %.0f%%' % (len(shell2), id2, len(shell2), 100*id2/len(shell2)))

# 5) superpose and per-residue side-chain RMSD
P = np.array([zA[h2z[k]]['ca'] for k in hk]); Q = np.array([hA[k]['ca'] for k in hk])
R, cP, cQ, rmsd = kabsch(P, Q)
print('global CA rmsd %.2f (%d pairs)' % (rmsd, len(hk)))
P3 = np.array([z3A[h2z[k]]['ca'] for k in hk]); R3, cP3, cQ3, rmsd3 = kabsch(P3, Q)
print('AF3 global CA rmsd %.2f' % rmsd3)

BB = {'N', 'CA', 'C', 'O', 'CB'}
def sc_rmsd(k):
    hr = hA[k]; zr = zA[h2z[k]]
    zmap = {}
    for a in zr['atoms']:
        zmap.setdefault(a['atom'], a)
    pairs = []
    for a in hr['atoms']:
        if a['atom'] in BB or a['atom'] not in zmap: continue
        pairs.append((zmap[a['atom']], a))
    if not pairs: return None
    Pm = np.array([[a['x'], a['y'], a['z']] for a, b in pairs])
    Qm = np.array([[b['x'], b['y'], b['z']] for a, b in pairs])
    Pm = (Pm - cP) @ R + cQ
    return math.sqrt(((Pm - Qm)**2).sum(1).mean())

rows = []
for k, dd in sorted(pocket):
    zk_ = h2z[k]
    haa = hA[k]['aa']; zaa = zA[zk_]['aa']
    rows.append(dict(pdb=k, unp=offsets[k], zf=zk_, h=haa, z=zaa,
                     ident=haa == zaa, dmin=dd,
                     plddt=zA[zk_]['b'], plddt3=z3A[zk_]['b'],
                     scrmsd=sc_rmsd(k), bfac=hA[k]['b']))
print('\n%4s %4s %4s %3s %3s %6s %6s %6s %6s %6s' % ('PDB', 'UNP', 'ZF', 'H', 'Z', 'dmin', 'pLDDT', 'pLDDT3', 'scRMSD', 'B'))
for r in rows:
    print('%4d %4d %4d %3s %3s %6.2f %6.1f %6.1f %6s %6.1f' % (
        r['pdb'], r['unp'], r['zf'], r['h'], r['z'], r['dmin'], r['plddt'], r['plddt3'],
        ('%.2f' % r['scrmsd']) if r['scrmsd'] is not None else '-', r['bfac']))

json.dump(dict(rows=rows), open(os.path.join(base, 'pocket_rows.json'), 'w'), indent=1)

# loop 470-480 backbone rmsd (zebrafish numbering)
lo = [k for k in hk if 470 <= h2z[k] <= 482]
Pm = np.array([(zA[h2z[k]]['ca'] - cP) @ R + cQ for k in lo])
Qm = np.array([hA[k]['ca'] for k in lo])
print('\nzebrafish 470-482 CA rmsd after global fit: %.2f A (%d residues)' % (np.sqrt(((Pm-Qm)**2).sum(1).mean()), len(lo)))

# terminal truncation info
print('\n9R6N chain A residues %d-%d; UniProt %d-%d' % (hk[0], hk[-1], offsets[hk[0]], offsets[hk[-1]]))
print('zebrafish model residues %d-%d' % (zk[0], zk[-1]))
