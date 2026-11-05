p="base_clean.pml"
s=open(p,encoding="utf-8").read()
n=s.count("0xCFEAF6")
s=s.replace("0xCFEAF6","0xF5C157")
open(p,"w",encoding="utf-8").write(s)
print("replaced",n)
