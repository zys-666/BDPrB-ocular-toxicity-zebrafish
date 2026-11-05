p="snap.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("zoom 0.95\npng snap/snap_0.png, dpi=200, ray=1","zoom 0.95\nviewport 800, 800\npng snap/snap_0.png, dpi=200, ray=1")
open(p,"w",encoding="utf-8").write(s)
print("patched")
