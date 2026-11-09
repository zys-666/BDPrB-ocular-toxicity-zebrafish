p="redock_overlay_toon.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("color 0x9FB0B8, pocket","color 0xA5B59A, pocket")
open(p,"w",encoding="utf-8").write(s)
print("patched")
