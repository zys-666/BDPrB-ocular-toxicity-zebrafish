import numpy as np, os
base = r"D:\user data\Desktop\斑马鱼\5.dockMD\MMPBSA"
d = np.loadtxt(os.path.join(base, "_pid~MMPBSA.dat"), comments="#")
t, dg = d[:,0], d[:,1]
print("n=%d  mean %.1f  median %.1f  SD %.1f" % (len(dg), dg.mean(), np.median(dg), dg.std()))
out = dg > -100
print("outliers (dG > -100 kJ/mol): %d frames" % out.sum())
print("   times:", np.round(t[out],1))
print("   values:", np.round(dg[out],1))
m = dg[~out]
print("excluding outliers: mean %.1f +/- %.1f (SEM), SD %.1f, n=%d" % (m.mean(), m.std(ddof=1)/len(m)**0.5, m.std(ddof=1), len(m)))
print("median-based 95%% CI: %.1f - %.1f" % (np.percentile(dg,2.5), np.percentile(dg,97.5)))
# components for an outlier frame
for i in np.where(out)[0][:3]:
    print("frame t=%.1f: dG=%.1f MM=%.1f COU=%.1f VDW=%.1f PB=%.1f SA=%.1f" % (t[i], dg[i], d[i,3], d[i,4], d[i,5], d[i,6], d[i,7]))
