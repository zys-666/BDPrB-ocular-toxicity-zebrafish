src="KyA33_crystal_vs_redock_in_protein.pdb"
import os
src=os.path.join(r"D:\user data\Desktop\斑马鱼\5.dockMD\dock\docking_results", src)
dst="KyA33_overlay_uniprot.pdb"
n=0
with open(dst,"w") as f:
    for l in open(src):
        if l.startswith(("ATOM","HETATM")) and l[21]=="A":
            rid=int(l[22:26])+23
            l=l[:22]+"%4d"%rid+l[26:]
            n+=1
        f.write(l)
print("renumbered", n, "protein atoms ->", dst)
