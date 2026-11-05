import numpy as np
def lig_and_sg_gro(path):
    L=open(path).read().splitlines(); n=int(L[1]); lig=[]; sg=None
    for l in L[2:2+n]:
        rn=l[5:10].strip(); rid=int(l[0:5]); an=l[10:15].strip()
        x=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
        if rn=='LIG': lig.append(x)
        if rid==320 and an=='SG' and rn=='CYS': sg=x
    return np.array(lig), sg
for p in ['../02_equil/em_whole.gro','../01_topology/complex.gro']:
    lig,sg=lig_and_sg_gro(p)
    print(p.split('/')[-1], 'lig centroid to C320 SG = %.2f A' % np.linalg.norm(lig.mean(0)-sg))
# docked pose
pose=[]
for l in open('../../dock/dock_run/dock_out.pdbqt'):
    if l.startswith('MODEL'): break
for l in open('../../dock/dock_run/dock_out.pdbqt'):
    if l.startswith(('ATOM','HETATM')): pose.append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    elif l.startswith('ENDMDL'): break
pose=np.array(pose)
sgd=None
for l in open('../../dock/dock_run/protein.pdb'):
    if l.startswith('ATOM') and int(l[22:26])==320 and l[12:16].strip()=='SG':
        sgd=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
print('docked pose centroid to C320 SG = %.2f A' % np.linalg.norm(pose.mean(0)-sgd))
