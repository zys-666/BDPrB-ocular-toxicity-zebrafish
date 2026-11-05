import numpy as np
from collections import OrderedDict
def read_gro_lig(path):
    L=open(path).read().splitlines(); n=int(L[1]); d={}
    for l in L[2:2+n]:
        if l[5:10].strip()=='LIG':
            d[l[10:15].strip()]=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
    return d
# protein CA for alignment
def cas(path):
    L=open(path).read().splitlines(); n=int(L[1]); d={}
    for l in L[2:2+n]:
        rn=l[5:10].strip(); rid=int(l[0:5]); an=l[10:15].strip()
        if an=='CA' and rn not in ('SOL','NA','CL','LIG'):
            d[rid]=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
    return d
E=cas('../02_equil/em_whole.gro'); D=cas('../../dock/dock_run/protein.pdb')
ks=sorted(set(E)&set(D))
P=np.array([E[k] for k in ks]); Q=np.array([D[k] for k in ks])
Pm,Qm=P.mean(0),Q.mean(0); Pc,Qc=P-Pm,Q-Qm
U,S,Vt=np.linalg.svd(Pc.T@Qc); dd=np.sign(np.linalg.det(U@Vt)); R=U@np.diag([1,1,dd])@Vt
lig=read_gro_lig('../02_equil/em_whole.gro')
X=np.array([lig[k] for k in lig]); Xd=(X-Pm)@R+Qm
print('EM ligand centroid in dock frame:', Xd.mean(0).round(2))
print('distance to KyA33 pocket centre (-12.34,3.55,-4.20): %.2f' % np.linalg.norm(Xd.mean(0)-np.array([-12.34,3.55,-4.20])))

