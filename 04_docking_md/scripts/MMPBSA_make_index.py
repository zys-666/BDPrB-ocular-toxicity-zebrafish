import os
src = r"D:\user data\Desktop\斑马鱼\5.dockMD\MD\03_production\index.ndx"
dst = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA\index.ndx"

groups = {}
order = []
name = None
for line in open(src):
    s = line.strip()
    if s.startswith('[') and s.endswith(']'):
        name = s[1:-1].strip()
        groups[name] = []
        order.append(name)
    elif name and s:
        groups[name].extend(int(x) for x in s.split())

print("groups:", order)
prot = groups["Protein"]
lig = groups["System"]  # placeholder
# find the ligand group (atoms of residue LIG) - the auto group is named LIG
lig = groups.get("LIG")
if lig is None:
    raise SystemExit("no LIG group")
combo = sorted(set(prot) | set(lig))

with open(dst, 'w') as f:
    def write(nm, atoms):
        f.write("[ %s ]\n" % nm)
        for i in range(0, len(atoms), 15):
            f.write(" ".join(str(a) for a in atoms[i:i+15]) + "\n")
    write("Protein", prot)
    write("LIGAND", lig)
    write("Protein_LIGAND", combo)
    write("System", groups["System"])
print("protein atoms:", len(prot), "ligand atoms:", len(lig), "complex:", len(combo))
