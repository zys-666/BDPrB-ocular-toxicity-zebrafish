import numpy as np
for f in ["test_a.pdb","test_b.pdb","test_c.pdb"]:
    prot=[]; lig=[]
    for l in open(f):
        if l.startswith(("ATOM","HETATM")):
            x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            (lig if l[17:20].strip()=="KYA" else prot).append(x)
    d=np.linalg.norm(np.array(lig)[:,None]-np.array(prot)[None],axis=-1).min()
    print("%-10s lig %d, min dist %.2f A" % (f,len(lig),d))
