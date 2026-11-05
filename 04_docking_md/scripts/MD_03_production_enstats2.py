import numpy as np
p=r"..\04_analysis\md_energy.xvg"
rows=[]
for l in open(p):
    if l.startswith(("#","@")) or not l.strip(): continue
    parts=l.split()
    if len(parts)==5:
        try: rows.append([float(x) for x in parts])
        except: pass
d=np.array(rows)
print("frames:", len(d), "time range %.0f-%.0f ps" % (d[:,0].min(), d[:,0].max()))
names=["Time","Temperature","Pressure","Density","Potential"]
for i in range(1,5):
    c=d[:,i]; m=d[:,0]>=d[:,0].max()*0.5
    print("  %-12s mean %.4g   last-50%% %.4g   std %.3g" % (names[i], c.mean(), c[m].mean(), c.std()))
