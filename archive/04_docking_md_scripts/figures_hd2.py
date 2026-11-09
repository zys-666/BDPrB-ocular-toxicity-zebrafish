p="base_clean.pml"
s=open(p,encoding="utf-8").read()
n=s.count("0xFBD9C9")
s=s.replace("0xFBD9C9","0xCFEAF6")
open(p,"w",encoding="utf-8").write(s)
print("replaced", n, "occurrences")
