# -*- coding: utf-8 -*-
import numpy as np, math, os, sys
from collections import OrderedDict

gro = sys.argv[1] if len(sys.argv) > 1 else "em.gro"
L = open(gro).read().splitlines()
n = int(L[1])
box = [float(x) for x in L[2+n].split()[:3]]

res = OrderedDict()
for l in L[2:2+n]:
    rn = l[5:10].strip(); rid = int(l[0:5]); an = l[10:15].strip()
    xyz = np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]) * 10.0
    res.setdefault((rid, rn), {})[an] = xyz

prot = [(k, v) for k, v in res.items() if k[1] not in ("SOL", "NA", "CL", "LIG")]
lig = [(k, v) for k, v in res.items() if k[1] == "LIG"]
print("atoms %d, residues %d (protein %d, ligand %d)" % (n, len(res), len(prot), len(lig)))

def dih(p1, p2, p3, p4):
    b0 = p1 - p2; b1 = p3 - p2; b2 = p4 - p3
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))

print("dihedral self-test (expect +90): %.1f" % dih(np.array([1.,0,0]), np.zeros(3), np.array([0.,1,0]), np.array([0.,1,1])))

cas = [v["CA"] for k, v in prot if "CA" in v]
ca_d = [np.linalg.norm(cas[i+1]-cas[i]) for i in range(len(cas)-1)]
print("CA-CA: mean %.2f, min %.2f, max %.2f, >4.5A: %d / %d" % (np.mean(ca_d), np.min(ca_d), np.max(ca_d), sum(d>4.5 for d in ca_d), len(ca_d)))

omegas = []
for i in range(len(prot)-1):
    a = prot[i][1]; b = prot[i+1][1]
    if all(k in a for k in ("CA","C")) and all(k in b for k in ("N","CA")):
        omegas.append(dih(a["CA"], a["C"], b["N"], b["CA"]))
omegas = np.array(omegas)
print("omega: mean|w| %.1f, min %.1f, max %.1f, cis(|w|<30): %d / %d" % (np.abs(omegas).mean(), omegas.min(), omegas.max(), (np.abs(omegas)<30).sum(), len(omegas)))

regions = dict(alpha=0, beta=0, left=0, other=0); nres = 0
for i in range(1, len(prot)-1):
    a, b, c = prot[i-1][1], prot[i][1], prot[i+1][1]
    if "C" not in a or not all(k in b for k in ("N","CA","C")) or "N" not in c: continue
    phi = dih(a["C"], b["N"], b["CA"], b["C"])
    psi = dih(b["N"], b["CA"], b["C"], c["N"])
    nres += 1
    if -160 <= phi <= -30 and -75 <= psi <= 30: regions["alpha"] += 1
    elif -180 <= phi <= -40 and (90 <= psi <= 180 or -180 <= psi <= -150): regions["beta"] += 1
    elif 20 <= phi <= 110 and -30 <= psi <= 110: regions["left"] += 1
    else: regions["other"] += 1
print("phi/psi: alpha %d (%.0f%%), beta %d (%.0f%%), left %d, other %d of %d" % (
    regions["alpha"], 100*regions["alpha"]/nres, regions["beta"], 100*regions["beta"]/nres,
    regions["left"], regions["other"], nres))

lig_atoms = lig[0][1]
print("ligand atoms:", len(lig_atoms))
itp = open(os.path.join("..", "01_topology", "LIG.itp")).read().splitlines()
i = itp.index("[ atoms ]"); names = []
for l in itp[i+1:]:
    if l.strip().startswith("["): break
    p = l.split()
    if len(p) >= 5 and not l.strip().startswith(";"): names.append(p[4])
i = itp.index("[ bonds ]"); bl = []
for l in itp[i+1:]:
    if l.strip().startswith("["): break
    p = l.split()
    if len(p) >= 4 and not l.strip().startswith(";"):
        bl.append(np.linalg.norm(lig_atoms[names[int(p[0])-1]] - lig_atoms[names[int(p[1])-1]]))
bl = np.array(bl)
print("ligand bonds: n=%d, mean %.3f A, min %.3f, max %.3f" % (len(bl), bl.mean(), bl.min(), bl.max()))

prot_xyz = np.array([x for k, v in prot for an, x in v.items() if not an.startswith("H")])
lig_xyz = np.array([x for an, x in lig_atoms.items() if not an.startswith("H")])
D = np.linalg.norm(lig_xyz[:, None, :] - prot_xyz[None, :, :], axis=-1)
print("ligand-protein heavy-atom min distance: %.2f A, contacts <4.0 A: %d" % (D.min(), (D < 4.0).sum()))
