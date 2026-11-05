import numpy as np
base=r"D:\user data\Desktop\斑马鱼\5.dockMD\dock\docking_results"
for f in ["complex_BDPrB.pdb","complex_KyA33.pdb"]:
    n=0; lig=0; chains=set()
    for l in open(base+"\\"+f):
        if l.startswith(("ATOM","HETATM")):
            n+=1; chains.add(l[21])
            if l[21]=="L": lig+=1
    print("%-22s %d atoms, ligand %d, chains %s" % (f,n,lig,sorted(chains)))
