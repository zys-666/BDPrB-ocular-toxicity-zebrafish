import glob, re
for p in ["base_clean.pml"] + glob.glob("probe_*.pml"):
    s=open(p,encoding="utf-8").read()
    s=s.replace("zoom (cplx and chain C), 11","zoom (cplx and chain C), 9")
    open(p,"w",encoding="utf-8").write(s)
print("zoom set to 9")
