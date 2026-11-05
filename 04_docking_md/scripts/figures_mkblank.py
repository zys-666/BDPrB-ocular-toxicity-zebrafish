s=open("hot_probe_195.pml",encoding="utf-8").read()
lines=[l for l in s.splitlines() if not l.startswith("label ")]
open("hot_blank.pml","w",encoding="utf-8").write("\n".join(lines)+"\n")
print("wrote hot_blank.pml")
