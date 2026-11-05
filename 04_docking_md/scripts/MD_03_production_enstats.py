import numpy as np, os
p=r"..\04_analysis\md_energy.xvg"
if os.path.exists(p):
    d=np.loadtxt(p, comments=("#","@"))
    print("columns:", d.shape[1]-1, "frames:", len(d))
    t=d[:,0]
    for i in range(1, d.shape[1]):
        c=d[:,i]
        print("  col%d: mean %.4g  last-50%% %.4g  std %.4g" % (i, c.mean(), c[t>=t.max()*0.5].mean(), c.std()))
else:
    print("missing")
