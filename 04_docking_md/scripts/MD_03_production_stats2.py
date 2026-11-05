import numpy as np, glob, os
for f in ['rmsd_bb.xvg','rmsd_lig.xvg','gyrate.xvg','rmsf.xvg']:
    p=os.path.join('..','04_analysis',f)
    if not os.path.exists(p):
        print(f,'MISSING'); continue
    d=np.loadtxt(p, comments=('#','@'))
    t=d[:,0]
    print('%s: %d frames, t=%.1f-%.1f ns' % (f, len(t), t.min(), t.max()))
    for i in range(1,d.shape[1]):
        c=d[:,i]
        half=c[t>=t.max()*0.5]
        print('   col%d: mean %.3f, last-50%% mean %.3f, max %.3f' % (i, c.mean(), half.mean(), c.max()))
