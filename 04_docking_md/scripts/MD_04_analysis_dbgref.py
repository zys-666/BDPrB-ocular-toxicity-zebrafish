import numpy as np, os
from collections import OrderedDict
lines=open('../../dock/dock_run/ligand.mol2').read().splitlines()
i=lines.index('@<TRIPOS>ATOM'); j=lines.index('@<TRIPOS>BOND')
atoms=[l.split() for l in lines[i+1:j]]
names=[a[1] for a in atoms]
bonds=[(int(b[1]),int(b[2])) for b in (l.split() for l in lines[j+1:] if l.strip() and not l.startswith('@'))]
models=[];cur=None
for l in open('../../dock/dock_run/dock_out.pdbqt'):
    if l.startswith('MODEL'): cur=[]
    elif l.startswith('ENDMDL'): models.append(np.array(cur)); cur=None
    elif l.startswith(('ATOM','HETATM')) and cur is not None:
        cur.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
pose=models[0]; perm=[0,1,10,12,13,15,16,17,4,7,5,3,6,8,9,11,14,2,18,19,20,21]
ref=np.zeros((22,3))
for k in range(22): ref[perm[k]]=pose[k]
bl=[np.linalg.norm(ref[a-1]-ref[b-1]) for a,b in bonds if a<=22 and b<=22]
print('ref heavy bond lengths: min %.2f max %.2f mean %.2f' % (min(bl),max(bl),np.mean(bl)))
# original mol2 coords for comparison
orig=np.array([[float(a[2]),float(a[3]),float(a[4])] for a in atoms[:22]])
print('rmsd ref vs original mol2 (heavy, no fit): %.2f' % np.sqrt(((ref-orig)**2).sum(1).mean()))
Qc=ref-ref.mean(0); Pc=orig-orig.mean(0)
U,S,Vt=np.linalg.svd(Pc.T@Qc); dd=np.sign(np.linalg.det(U@Vt)); R=U@np.diag([1,1,dd])@Vt
print('rmsd ref vs original mol2 (fit): %.2f' % np.sqrt((((Pc@R)-Qc)**2).sum(1).mean()))
# internal torsions compare: list 7 rotatable bonds
print('ref centroid', ref.mean(0).round(2))
