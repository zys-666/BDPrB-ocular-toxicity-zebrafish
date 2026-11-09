import numpy as np
def ca(path):
    d={}
    for l in open(path):
        if l.startswith(('ATOM','HETATM')) and l[12:16].strip()=='CA':
            d[int(l[22:26])]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    return d
for f in ['chk50.pdb','chk100.pdb']:
    d=ca(f); ks=sorted(d)
    dist=[np.linalg.norm(d[ks[i+1]]-d[ks[i]]) for i in range(len(ks)-1)]
    print('%s: %d CA, CA-CA mean %.2f min %.2f max %.2f, >4.5A: %d' % (
        f, len(ks), np.mean(dist), np.min(dist), np.max(dist), sum(x>4.5 for x in dist)))
