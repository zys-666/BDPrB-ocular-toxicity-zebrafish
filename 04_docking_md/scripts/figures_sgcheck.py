import numpy as np
path="../dock/docking_results/complex_BDPrB.pdb"
SG=None; L=[]
for l in open(path):
    if not l.startswith(("ATOM","HETATM")): continue
    x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    if l[21]=="L":
        L.append((l[12:16].strip(), x))
    elif l[17:20].strip()=="CYS" and l[22:26].strip()=="320" and l[12:16].strip()=="SG":
        SG=x
for name,x in sorted(L, key=lambda t: np.linalg.norm(t[1]-SG))[:5]:
    print("%-4s  %.2f A" % (name, np.linalg.norm(x-SG)))
