import numpy as np, glob
for f in ['rmsd_backbone.xvg','rmsd_ligand.xvg','gyrate.xvg']:
    try:
        d=np.loadtxt(f, comments=('#','@'))
        print(f)
        for i in range(1,d.shape[1]):
            col=d[:,i]
            print('   col%d: mean %.3f  last %.3f  min %.3f max %.3f' % (i, col.mean(), col[-1], col.min(), col.max()))
    except Exception as e:
        print(f,'ERR',e)
