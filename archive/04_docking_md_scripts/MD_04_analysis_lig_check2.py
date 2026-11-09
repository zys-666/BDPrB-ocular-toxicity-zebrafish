import numpy as np, os
from collections import OrderedDict
def read_gro(path):
    L=open(path).read().splitlines(); n=int(L[1]); d=OrderedDict()
    for l in L[2:2+n]:
        rn=l[5:10].strip(); rid=int(l[0:5]); an=l[10:15].strip()
        xyz=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
        d.setdefault((rid,rn),{})[an]=xyz
    return d
d=read_gro('../02_equil/em_whole.gro')
lig=[v for k,v in d.items() if k[1]=='LIG'][0]
prot_ca={rid:at['CA'] for (rid,rn),at in d.items() if rn not in ('SOL','NA','CL','LIG') and 'CA' in at}
# receptor reference (docked frame)
pdb={}
for l in open('../../dock/dock_run/protein.pdb'):
    if l.startswith('ATOM') and l[12:16].strip()=='CA':
        pdb[int(l[22:26])]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
ks=sorted(set(prot_ca)&set(pdb))
P=np.array([prot_ca[k] for k in ks]); Q=np.array([pdb[k] for k in ks])
Pc=P-P.mean(0); Qc=Q-Q.mean(0)
U,S,Vt=np.linalg.svd(Pc.T@Qc); dd=np.sign(np.linalg.det(U@Vt)); R=U@np.diag([1,1,dd])@Vt
rms=np.sqrt(((Pc@R-Qc)**2).sum()/len(P))
print('protein CA fit RMSD em vs dock: %.2f A (%d CA)' % (rms,len(P)))
t=Q.mean(0)-P.mean(0)@R
lig_in_dock=(np.array([lig[n] for n in lig])-P.mean(0))@R+Q.mean(0)
# reference pose
lines=open('../../dock/dock_run/ligand.mol2').read().splitlines()
i=lines.index('@<TRIPOS>ATOM'); j=lines.index('@<TRIPOS>BOND')
names=[l.split()[1] for l in lines[i+1:j]]
models=[];cur=None
for l in open('../../dock/dock_run/dock_out.pdbqt'):
    if l.startswith('MODEL'): cur=[]
    elif l.startswith('ENDMDL'): models.append(np.array(cur)); cur=None
    elif l.startswith(('ATOM','HETATM')) and cur is not None:
        cur.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
pose=models[0]; perm=[0,1,10,12,13,15,16,17,4,7,5,3,6,8,9,11,14,2,18,19,20,21]
ref=np.zeros((22,3))
for k in range(22): ref[perm[k]]=pose[k]
eml=np.array([lig[n] for n in names[:22]])
eml=(eml-P.mean(0))@R+Q.mean(0)
print('ligand heavy-atom RMSD after protein-based alignment: %.2f A' % np.sqrt(((eml-ref)**2).sum(1).mean()))
print('ligand centroid shift: %.2f A' % np.linalg.norm(eml.mean(0)-ref.mean(0)))
for n,c in [('CYS320','SG'),('TRP195','CZ2'),('MET192','SD')]:
    rid=320 if n=='CYS320' else (195 if n=='TRP195' else 192)
    x=d[(rid,n)]['SG' if n=='CYS320' else ('CZ2' if n=='TRP195' else 'SD')]
    x=(x-P.mean(0))@R+Q.mean(0)
    print('ligand-%s:%s min dist %.2f A' % (n,c, np.min(np.linalg.norm(ref-x,axis=1))))
