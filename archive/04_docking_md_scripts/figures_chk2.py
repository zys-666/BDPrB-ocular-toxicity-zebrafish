import numpy as np
for t in [0,25000,50000,75000,100000]:
    P=[];L=[]
    for l in open("snap/frame_%d.pdb"%t):
        if l.startswith(("ATOM","HETATM")):
            x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            (L if l[21]=="B" else P).append(x)
    P=np.array(P);L=np.array(L)
    print("t=%6d ps  ligand %d atoms, centroid dist to protein %.1f A, min dist %.2f A" % (
        t, len(L), np.linalg.norm(L.mean(0)-P.mean(0)), np.linalg.norm(L[:,None]-P[None],axis=-1).min()))
