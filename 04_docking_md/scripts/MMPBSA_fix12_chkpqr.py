import numpy as np, glob, os
base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA\fix12"
for f in glob.glob(os.path.join(base,'fix*','*_com.pqr'))[:3]:
    lig=[]; prot=[]
    for l in open(f):
        if l.startswith(('ATOM','HETATM')):
            x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            (lig if l[17:20].strip() in ('LIG','UNL') else prot).append(x)
    if lig and prot:
        d=np.linalg.norm(np.array(lig)[:,None]-np.array(prot)[None],axis=-1).min()
        print(os.path.basename(f), 'ligand %d atoms, min dist %.2f A' % (len(lig), d))
