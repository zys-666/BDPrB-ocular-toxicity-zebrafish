p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\redock_overlay_toon.pml"
s=open(p,encoding='utf-8').read()
s=s.replace("set stick_radius, 0.13, pocket",
"""set stick_radius, 0.13, pocket
set cartoon_color, 0xFBD9C9, cplx and chain A""")
open(p,'w',encoding='utf-8').write(s)
print("patched")
