import numpy as np
base = r"D:\user data\Desktop\斑马鱼\5.dockMD\controls\apo"
for f in ["rmsd_bb.xvg","gyrate.xvg"]:
    d=np.loadtxt(base+"\\"+f, comments=("#","@"))
    t=d[:,0]; c=d[:,1]
    m=t>=50
    print("%s: mean %.3f nm, last-50ns %.3f nm, max %.3f" % (f, c.mean(), c[m].mean(), c.max()))
