import numpy as np
p=r"D:\user data\Desktop\斑马鱼\5.dockMD\dock\docking_results\KyA33_crystal_vs_redock_in_protein.pdb"
prot=[];C=[];D=[]
for l in open(p):
    if l.startswith(("ATOM","HETATM")):
        x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
        ch=l[21]
        if ch=="A": prot.append(x)
        elif ch=="C": C.append(x)
        elif ch=="D": D.append(x)
P=np.array(prot);Cc=np.array(C);Dd=np.array(D)
print("protein %d, crystal %d, redock %d" % (len(P),len(Cc),len(Dd)))
print("min dist crystal-protein %.2f A, redock-protein %.2f A" % (
    np.linalg.norm(Cc[:,None]-P[None],axis=-1).min(), np.linalg.norm(Dd[:,None]-P[None],axis=-1).min()))
print("centroid crystal-redock %.2f A" % np.linalg.norm(Cc.mean(0)-Dd.mean(0)))
