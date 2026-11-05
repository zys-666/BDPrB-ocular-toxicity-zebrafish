import os
base=r"D:\user data\Desktop\斑马鱼\5.dockMD"
cry=os.path.join(base,"dock","recontrol","lig_ref.pdb")
red=os.path.join(base,"dock","docking_results","KyA33_redock_pose1.pdb")
out=os.path.join(base,"dock","docking_results","KyA33_crystal_vs_redock.pdb")
n=0
with open(out,"w") as f:
    for l in open(cry):
        if l.startswith(("ATOM","HETATM")):
            n+=1
            f.write(l[:6]+"%5d"%n+l[11:21]+"C"+l[22:])
    f.write("TER\n")
    for l in open(red):
        if l.startswith(("ATOM","HETATM")):
            n+=1
            f.write(l[:6]+"%5d"%n+l[11:21]+"D"+l[22:])
    f.write("END\n")
print("wrote",out,n,"atoms")
