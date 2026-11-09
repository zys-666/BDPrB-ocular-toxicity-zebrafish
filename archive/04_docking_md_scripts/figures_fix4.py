p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\redock_overlay_bead.pml"
s=open(p,encoding='utf-8').read()
s=s.replace("set transparency, 0.55, prot",
"""# outline ray mode does not support object transparency -> cut a window instead
select near, byres (prot within 9.0 of ligC)
hide spheres, near""")
open(p,'w',encoding='utf-8').write(s)
print("patched")
