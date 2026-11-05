import numpy as np
L=open('../01_topology/complex.gro').read().splitlines(); n=int(L[1])
prot=[]; lig=[]
for l in L[2:2+n]:
    rn=l[5:10].strip(); an=l[10:15].strip()
    x=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
    if rn=='LIG': lig.append((an,x))
    elif rn not in ('SOL','NA','CL'): prot.append((an,x))
Lg=np.array([x for a,x in lig if not a.startswith('H')])
P=np.array([x for a,x in prot if not a.startswith('H')])
D=np.linalg.norm(Lg[:,None]-P[None],axis=-1)
print('complex.gro: ligand-protein heavy min dist %.2f A, contacts<3.0: %d' % (D.min(), (D<3.0).sum()))
# any protein atom (incl H) closer than 2.0?
Pall=np.array([x for a,x in prot]); Dall=np.linalg.norm(Lg[:,None]-Pall[None],axis=-1)
print('incl. H: min %.2f A' % Dall.min())
