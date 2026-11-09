p="snap5.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("set auto_zoom, 0\ndisable f1 f2 f3 f4",
"""python
from pymol import cmd
for o in ["f1", "f2", "f3", "f4"]:
    cmd.align(o + " and name CA", "f0 and name CA")
python end

set auto_zoom, 0
disable f1 f2 f3 f4""")
open(p,"w",encoding="utf-8").write(s)
print("patched")
