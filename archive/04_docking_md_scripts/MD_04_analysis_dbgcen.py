import numpy as np
L=open('../02_equil/em_whole.gro').read().splitlines(); n=int(L[1])
prot=[]; lig=[]
for l in L[2:2+n]:
    rn=l[5:10].strip()
    x=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
    if rn=='LIG': lig.append(x)
    elif rn not in ('SOL','NA','CL'): prot.append(x)
prot=np.array(prot); lig=np.array(lig)
print('protein centroid', prot.mean(0).round(1), 'ligand centroid', lig.mean(0).round(1))
D=np.linalg.norm(lig[:,None]-prot[None],axis=-1)
print('raw min dist %.2f' % D.min())
box=np.array([float(x) for x in L[2+n].split()])
print('box', box)
