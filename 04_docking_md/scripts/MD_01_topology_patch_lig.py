p='LIG.itp'
s=open(p).read()
if 'POSRES_LIG' not in s:
    s=s.rstrip()+'\n\n; ligand position restraints\n#ifdef POSRES_LIG\n#include "ligand/posre_LIG.itp"\n#endif\n'
    open(p,'w').write(s)
print('ligand itp patched')
