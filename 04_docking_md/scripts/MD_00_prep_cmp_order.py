pqt=[l for l in open(r'D:\user data\Desktop\斑马鱼\5.dockMD\dock\dock_run\ligand.pdbqt') if l.startswith(('ATOM','HETATM'))]
pe=[l[77:79].strip() for l in pqt]
print('PDBQT elements:', pe)
lines=open(r'D:\user data\Desktop\斑马鱼\5.dockMD\dock\dock_run\ligand.mol2').read().splitlines()
i=lines.index('@<TRIPOS>ATOM'); j=lines.index('@<TRIPOS>BOND')
me=[]
for l in lines[i+1:j]:
    p=l.split()
    me.append(p[5].split('.')[0])
print('MOL2  elements:', me[:22])
print('counts pdbqt:', {e:pe.count(e) for e in set(pe)})
print('counts mol2 :', {e:me[:22].count(e) for e in set(me[:22])})
