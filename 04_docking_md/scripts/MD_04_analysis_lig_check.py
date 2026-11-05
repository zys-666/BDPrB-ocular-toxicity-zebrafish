import numpy as np
from collections import OrderedDict
# ---- ligand pose in em_whole vs docked pose ----
def read_gro(path):
    L=open(path).read().splitlines(); n=int(L[1]); d=OrderedDict()
    for l in L[2:2+n]:
        rn=l[5:10].strip(); rid=int(l[0:5]); an=l[10:15].strip()
        xyz=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
        d.setdefault((rid,rn),{})[an]=xyz
    return d
d=read_gro('../02_equil/em_whole.gro')
lig=[v for k,v in d.items() if k[1]=='LIG'][0]
# docked pose (PDBQT model1 with mol2 order)
import os
dock='../../dock/dock_run'
lines=open(os.path.join(dock,'ligand.mol2')).read().splitlines()
i=lines.index('@<TRIPOS>ATOM'); j=lines.index('@<TRIPOS>BOND')
names=[l.split()[1] for l in lines[i+1:j]]
models=[];cur=None
for l in open(os.path.join(dock,'dock_out.pdbqt')):
    if l.startswith('MODEL'): cur=[]
    elif l.startswith('ENDMDL'): models.append(np.array(cur)); cur=None
    elif l.startswith(('ATOM','HETATM')) and cur is not None:
        cur.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
pose=models[0]
# permutation pdbqt->mol2 computed earlier
perm=[0,1,10,12,13,15,16,17,4,7,5,3,6,8,9,11,14,2,18,19,20,21]
ref=np.zeros((22,3))
for k in range(22): ref[perm[k]]=pose[k]
em=np.array([lig[names[a]] for a in range(22)])
print('ligand heavy-atom RMSD (no fit) vs docked pose: %.2f A' % np.sqrt(((em-ref)**2).sum(1).mean()))
emc=em-em.mean(0); rc=ref-ref.mean(0)
U,S,Vt=np.linalg.svd(emc.T@rc); dd=np.sign(np.linalg.det(U@Vt)); R=U@np.diag([1,1,dd])@Vt
print('ligand best-fit RMSD: %.2f A' % np.sqrt((((emc@R)-rc)**2).sum(1).mean()))
print('centroid shift: %.2f A' % np.linalg.norm(em.mean(0)-ref.mean(0)))
# ---- pocket residues ----
prot=[(k,v) for k,v in d.items() if k[1] not in ('SOL','NA','CL','LIG')]
res={}
for (rid,rn),at in prot:
    for an,x in at.items():
        if an.startswith('H'): continue
        dmin=min(np.linalg.norm(x-lig[n]) for n in ['C','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','F','F1','O'])
        if dmin<4.0: res[(rid,rn)]=min(res.get((rid,rn),9),dmin)
print('pocket residues (<4.0 A):', len(res))
print('  ', ', '.join('%s%d'%(rn,rid) for (rid,rn) in sorted(res)))
