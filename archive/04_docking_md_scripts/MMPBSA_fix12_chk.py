import numpy as np
def split(path):
    prot=[]; lig=[]
    for l in open(path):
        if l.startswith(('ATOM','HETATM')):
            x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            (lig if l[17:20].strip()=='LIG' else prot).append(x)
    return np.array(prot), np.array(lig)
for f in ['t_plain.pdb','t_pbc.pdb']:
    try:
        p,l=split(f)
        d=np.linalg.norm(l[:,None]-p[None],axis=-1).min()
        print('%-12s protein %d atoms, ligand %d, min dist %.2f A' % (f,len(p),len(l),d))
    except Exception as e:
        print(f,'ERR',e)
