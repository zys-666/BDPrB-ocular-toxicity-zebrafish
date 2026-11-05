import numpy as np
lines=open('lig_mapped.mol2').read().splitlines()
i=lines.index('@<TRIPOS>ATOM'); j=lines.index('@<TRIPOS>BOND')
at=[(l.split()[1], np.array([float(l.split()[2]),float(l.split()[3]),float(l.split()[4])]), l.split()[5]) for l in lines[i+1:j] if l.strip()]
X=np.array([a[1] for a in at])
D=np.linalg.norm(X[:,None]-X[None,:],axis=-1); np.fill_diagonal(D,99)
print('atoms',len(at),'min dist %.3f'%D.min())
bad=np.argwhere(D<0.9)
print('close pairs:',len(bad)//2)
for a,b in bad[:5]:
    if a<b: print('  ',at[a][0],at[b][0],round(D[a,b],2))
