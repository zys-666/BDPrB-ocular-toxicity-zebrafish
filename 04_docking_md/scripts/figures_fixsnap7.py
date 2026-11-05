p="snap5.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("disable f1 f2 f3 f4\norient f0","set auto_zoom, 0\ndisable f1 f2 f3 f4\norient f0")
open(p,"w",encoding="utf-8").write(s)
print("patched")
