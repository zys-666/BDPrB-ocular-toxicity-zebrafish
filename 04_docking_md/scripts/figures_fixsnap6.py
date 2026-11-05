p="snap5.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("center (f0 and chain B)\nzoom (f0 and chain B), 13",
            "orient f0\nzoom (f0 and chain B), 11")
open(p,"w",encoding="utf-8").write(s)
print("patched")
