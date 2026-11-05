import numpy as np, os
base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA\fix12"
def split(path):
    prot=[]; lig=[]
    for l in open(path):
        if l.startswith(('ATOM','HETATM')):
            x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            (lig if l[17:20].strip()=='LIG' else prot).append(x)
    return np.array(prot), np.array(lig)
for f in ['t_plain.pdb','t_pbc.pdb']:
    p=os.path.join(base,f)
    if not os.path.exists(p):
        print(f,'missing'); continue
    pr,li=split(p)
    d=np.linalg.norm(li[:,None]-pr[None],axis=-1).min()
    print('%-12s protein %4d, ligand %2d, min dist %6.2f A' % (f,len(pr),len(li),d))
