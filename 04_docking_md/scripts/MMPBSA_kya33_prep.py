# -*- coding: utf-8 -*-
"""Prepare the MM-PBSA run for the KyA33 control system."""
import os, shutil

root = r"D:\user data\Desktop\斑马鱼\5.dockMD"
src_script = os.path.join(root, "mmpbsa", "gmx_mmpbsa_WSL.bash")
if not os.path.exists(src_script):
    src_script = os.path.join(r"D:\Gromacs25_1_win\mmpbsa\gmxmmpbsa", "gmx_mmpbsa_WSL.bash")
dst = os.path.join(root, "MMPBSA_kya33")
os.makedirs(dst, exist_ok=True)

shutil.copy(src_script, os.path.join(dst, "gmx_mmpbsa.bash"))
for f in ["run_mmpbsa_parallel.sh", "merge_parallel.sh", "status.sh"]:
    shutil.copy(os.path.join(r"D:\Gromacs25_1_win\mmpbsa\gmxmmpbsa", f), dst)
shutil.copy(os.path.join(root, "controls", "kya33", "md.tpr"), dst)
shutil.copy(os.path.join(root, "controls", "kya33", "index.ndx"), dst)
shutil.copy(os.path.join(root, "controls", "kya33", "trj_50-100ns.xtc"), dst)

p = os.path.join(dst, "gmx_mmpbsa.bash")
s = open(p, encoding="utf-8", errors="replace").read()
s = s.replace("tpr=md_0_1.tpr", "tpr=md.tpr")
s = s.replace("ndx=index.ndx", "ndx=index.ndx")
s = s.replace("pro=Protein", "pro=Protein")
s = s.replace("lig=LIGAND", "lig=LIGAND")
s = s.replace('trjconv="$gmx trjconv -dt 300"', 'trjconv="$gmx trjconv"')
old = """	echo -e "$lig\\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $trjcnt &>>$err -pbc mol -center
	echo -e "$com\\n$com" | $trjconv  -s $tpx -n $idx -f $trjcnt -o $trjcls &>>$err -pbc cluster
	echo -e "$lig\\n$com" | $trjconv  -s $tpx -n $idx -f $trjcls -o $pdb    &>>$err -fit rot+trans"""
new = """	cp $trjwho $trjcnt
	cp $trjcnt $trjcls
	echo "$com" | $trjconv -s $tpx -n $idx -f $trjcls -o $pdb &>>$err"""
n = s.count(old)
s = s.replace(old, new)
open(p, "w", encoding="utf-8", newline="\n").write(s)
print("script patched, preprocessing replacements:", n)
print("files:", sorted(os.listdir(dst)))
