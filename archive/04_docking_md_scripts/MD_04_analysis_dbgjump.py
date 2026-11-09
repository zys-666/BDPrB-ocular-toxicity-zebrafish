import numpy as np
from collections import OrderedDict
L=open("../02_equil/em.gro").read().splitlines(); n=int(L[1])
res=OrderedDict()
for l in L[2:2+n]:
    rn=l[5:10].strip(); rid=int(l[0:5]); an=l[10:15].strip()
    xyz=np.array([float(l[20:28]),float(l[28:36]),float(l[36:44])])*10.0
    res.setdefault((rid,rn),{})[an]=xyz
prot=[(k,v) for k,v in res.items() if k[1] not in ("SOL","NA","CL","LIG")]
cas=[(k, v["CA"]) for k,v in prot if "CA" in v]
for i in range(len(cas)-1):
    d=np.linalg.norm(cas[i+1][1]-cas[i][1])
    if d>4.5: print("jump %.1f A between %s and %s" % (d, cas[i][0], cas[i+1][0]))
print("first", cas[0][0], "last", cas[-1][0], "count", len(cas))
print("residue ids continuous:", [cas[i][0][0] for i in range(len(cas))][:5], "...", [c[0][0] for c in cas][-5:])
