import numpy as np
base = r"D:\user data\Desktop\斑马鱼\5.dockMD\figures"
for f in ["ill_zfish_vs_human_superposed.pdb","ill_zfish_vs_human_sidebyside.pdb"]:
    A=[];B=[]
    for l in open(base+"\\"+f):
        if l.startswith("ATOM"):
            x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            (A if l[21]=="A" else B).append(x)
    A=np.array(A);B=np.array(B)
    print("%-40s A %d atoms, B %d atoms, centroid distance %.1f A" % (f,len(A),len(B),np.linalg.norm(A.mean(0)-B.mean(0))))
