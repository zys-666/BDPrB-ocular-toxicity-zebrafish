import os
p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\redock_overlay.pml"
s=open(p,encoding='utf-8').read()
s=s.replace("""pseudoatom lbl, pos=[0, 0, 0]
hide everything, lbl
label lbl, "redock RMSD = 0.52 A"
set label_size, 26
set label_color, black
set label_position, [0, 0, 0]
label lbl""",
"""pseudoatom lbl, selection=crystal
hide everything, lbl
label lbl, "redock RMSD = 0.52 A"
set label_size, 26
set label_color, black
set label_position, [0, 6, 0]
label lbl""")
open(p,'w',encoding='utf-8').write(s)
print("patched")
