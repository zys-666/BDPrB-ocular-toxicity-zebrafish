import numpy as np
d={}; allx=[]
for l in open('fit_chk.pdb'):
    if l.startswith(('ATOM','HETATM')):
        x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]); allx.append(x)
        if l[12:16].strip()=='CA': d[int(l[22:26])]=x
ks=sorted(d); dist=[np.linalg.norm(d[ks[i+1]]-d[ks[i]]) for i in range(len(ks)-1)]
X=np.array(allx)
print('frame 75 ns: protein centroid', X.mean(0).round(1), ' CA-CA max %.2f, breaks %d' % (max(dist), sum(x>4.5 for x in dist)))
