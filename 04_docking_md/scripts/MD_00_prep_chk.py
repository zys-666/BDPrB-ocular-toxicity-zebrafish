import numpy as np
lines=open('lig_pose.mol2').read().splitlines()
i=lines.index('@<TRIPOS>ATOM'); j=lines.index('@<TRIPOS>BOND')
at=[]
for l in lines[i+1:j]:
    p=l.split()
    if len(p)>=6:
        at.append((p[1], np.array([float(p[2]),float(p[3]),float(p[4])]), p[5]))
print('atoms', len(at))
X=np.array([a[1] for a in at])
D=np.linalg.norm(X[:,None,:]-X[None,:,:],axis=-1)
np.fill_diagonal(D, 99)
print('min pair distance:', D.min())
idx=np.argwhere((D<0.9))
print('close pairs (<0.9 A):', len(idx)//2)
for a,b in idx[:10]:
    if a<b: print('  ', at[a][0], at[b][0], round(D[a,b],3))
print('coordinate range:', X.min(0).round(2), X.max(0).round(2))
