import numpy as np
def load(path, lig_chain):
    P=[];L=[];SG=None
    for l in open(path):
        if not l.startswith(("ATOM","HETATM")): continue
        x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        if l[21]==lig_chain: L.append(x)
        else:
            P.append(x)
            if l[17:20].strip()=="CYS" and l[22:26].strip()=="320" and l[12:16].strip()=="SG":
                SG=x
    return np.array(P),np.array(L),SG
for tag,path,ch in [("BDPrB","../dock/docking_results/complex_BDPrB.pdb","L"),
                    ("KyA33(zf)","../dock/docking_results/complex_KyA33.pdb","L")]:
    P,L,SG=load(path,ch)
    d=np.linalg.norm(L[:,None]-P[None],axis=-1)
    print("%-10s lig atoms=%d  min dist to protein=%.2f A   min dist to Cys320 SG=%s" % (
        tag,len(L),d.min(), "%.2f A"%np.linalg.norm(L-SG,axis=1).min() if SG is not None else "NA"))
