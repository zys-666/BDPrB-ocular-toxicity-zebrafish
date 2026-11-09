import numpy as np
for t in [0,25000,50000,75000,100000]:
    p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\snap\frame_%d.pdb"%t
    P=[];L=[]
    for l in open(p):
        if l.startswith(("ATOM","HETATM")):
            x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            (L if l[21]=="L" else P).append(x)
    P=np.array(P);L=np.array(L)
    print("t=%6d ns  min dist %.2f A" % (t/1000, np.linalg.norm(L[:,None]-P[None],axis=-1).min()))
