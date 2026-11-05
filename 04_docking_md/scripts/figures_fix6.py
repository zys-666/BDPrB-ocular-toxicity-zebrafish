p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\redock_overlay_toon.pml"
s=open(p,encoding='utf-8').read()
s=s.replace("color 0xE3BCAE, pocket","color 0x9FB0B8, pocket")
s=s.replace("set stick_radius, 0.10, pocket","set stick_radius, 0.13, pocket")
open(p,'w',encoding='utf-8').write(s)
print("patched")
