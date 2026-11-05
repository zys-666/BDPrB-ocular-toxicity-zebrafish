import numpy as np, glob, os
base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA\fix12"
for f in sorted(glob.glob(os.path.join(base,'fix*','*_com.pqr')))[:4]:
    lig=[]; prot=[]
    for l in open(f):
        if l.startswith(('ATOM','HETATM')):
            p=l.split()
            if len(p) < 9: continue
            x=np.array([float(p[5]),float(p[6]),float(p[7])])
            (lig if p[3] in ('LIG','UNL') else prot).append(x)
    if lig and prot:
        d=np.linalg.norm(np.array(lig)[:,None]-np.array(prot)[None],axis=-1).min()
        print(os.path.basename(f), 'lig %d, prot %d, min dist %.2f A' % (len(lig), len(prot), d))
