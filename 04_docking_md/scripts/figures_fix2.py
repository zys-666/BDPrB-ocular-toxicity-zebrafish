p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\redock_in_protein.pml"
s=open(p,encoding='utf-8').read()
s=s.replace("""select pocket, byres ((cplx and chain A) within 4.0 of (cplx and chain C))
show lines, pocket
color grey55, pocket
set line_width, 1.0""",
"""select pocket, byres ((cplx and chain A) within 4.0 of (cplx and chain C))
show sticks, pocket
color grey70, pocket
set stick_radius, 0.08, pocket
hide (hydro and pocket)""")
open(p,'w',encoding='utf-8').write(s)
print("patched")
