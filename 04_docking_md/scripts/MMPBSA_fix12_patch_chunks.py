# -*- coding: utf-8 -*-
"""Patch the per-chunk gmx_mmpbsa.bash so it does NOT re-run its own PBC clustering
(the input frames are already whole and in the same image)."""
import os, glob

base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA\fix12"
old = """	echo -e "$lig\\n$com" | $trjconv  -s $tpx -n $idx -f $trjwho -o $trjcnt &>>$err -pbc mol -center
	echo -e "$com\\n$com" | $trjconv  -s $tpx -n $idx -f $trjcnt -o $trjcls &>>$err -pbc cluster
	echo -e "$lig\\n$com" | $trjconv  -s $tpx -n $idx -f $trjcls -o $pdb    &>>$err -fit rot+trans"""

new = """	cp $trjwho $trjcnt
	cp $trjcnt $trjcls
	echo "$com" | $trjconv -s $tpx -n $idx -f $trjcls -o $pdb &>>$err"""

for f in glob.glob(os.path.join(base, "fix?", "gmx_mmpbsa.bash")):
    s = open(f, encoding="utf-8", errors="replace").read()
    if old in s:
        s = s.replace(old, new)
        open(f, "w", encoding="utf-8", newline="\n").write(s)
        print("patched", f)
    else:
        print("PATTERN NOT FOUND in", f)
