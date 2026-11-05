import numpy as np
d=np.loadtxt("mindist.xvg", comments=("#","@"))
t=d[:,0]/1000.0; m=d[:,1]*10
print("ligand-protein min dist: mean %.2f A, max %.2f A" % (m.mean(), m.max()))
for thr in (5,10,20,30):
    print("  frames with min dist > %2d A: %d (%.0f%%)" % (thr,(m>thr).sum(),100*(m>thr).mean()))
i=int(np.argmax(m)); print("max at %.1f ns: %.1f A" % (t[i], m[i]))
