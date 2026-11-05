L=open("../01_topology/complex.gro").read().splitlines(); n=int(L[1])
for l in L[2:2+n]:
    rid=int(l[0:5]); an=l[10:15].strip()
    if rid in (159,160,161,162,509,510) and an=="CA":
        print(repr(l))
