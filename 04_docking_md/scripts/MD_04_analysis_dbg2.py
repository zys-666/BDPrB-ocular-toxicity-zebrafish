L=open("../02_equil/em.gro").read().splitlines(); n=int(L[1])
for l in L[2:2+n]:
    rid=int(l[0:5]); rn=l[5:10].strip(); an=l[10:15].strip()
    if rid in (159,160,161,162,509,510) and an=="CA":
        print(repr(l))
