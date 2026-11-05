import numpy as np
from collections import OrderedDict
def read_gro(path):
    L=open(path).read().splitlines(); n=int(L[1]); d=[]
    for l in L[2:2+n]:
        rn=l[5:10].strip(); rid=int(l[0:5]); an=l[10:15].strip()
        xyz=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
        if rn=='LIG': d.append((rid,an,xyz))
    return d
c=read_gro('../01_topology/complex.gro')
e=read_gro('../02_equil/em_whole.gro')
C=np.array([x for _,_,x in c]); E=np.array([x for _,_,x in e])
print('complex.gro ligand centroid', C.mean(0).round(2), 'n', len(C))
print('em_whole   ligand centroid', E.mean(0).round(2), 'n', len(E))
print('shift before EM (complex vs dock ref -13.3 3.0 -3.05): %.2f' % np.linalg.norm(C.mean(0)-np.array([-13.3,3.01,-3.05])))
print('shift during EM: %.2f A' % np.linalg.norm(E.mean(0)-C.mean(0)))
print('pairwise internal rmsd complex vs em (fit):')
Pc=C-C.mean(0);Qc=E-E.mean(0)
U,S,Vt=np.linalg.svd(Pc.T@Qc);dd=np.sign(np.linalg.det(U@Vt));R=U@np.diag([1,1,dd])@Vt
print('  %.2f A' % np.sqrt((((Pc@R)-Qc)**2).sum(1).mean()))
