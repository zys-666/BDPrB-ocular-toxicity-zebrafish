p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\redock_overlay_toon.pml"
s=open(p,encoding='utf-8').read()
s=s.replace("color palegreen, cplx and chain A","color 0xFBD9C9, cplx and chain A")
s=s.replace("color wheat, pocket","color 0xE3BCAE, pocket")
open(p,'w',encoding='utf-8').write(s)
print("patched")
