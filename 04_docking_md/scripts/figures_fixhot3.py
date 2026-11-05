p="gen_hot.py"
s=open(p,encoding="utf-8").read()
s=s.replace("zoom (cplx and chain L), 11","zoom (cplx and chain L), 7.5")
open(p,"w",encoding="utf-8").write(s)
print("patched")
