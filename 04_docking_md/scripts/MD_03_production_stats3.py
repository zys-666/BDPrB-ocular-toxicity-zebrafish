import numpy as np
def load(p):
    return np.loadtxt(p, comments=('#','@'))
d=load('mindist.xvg'); t=d[:,0]
print('mindist LIG-Protein: mean %.2f nm, min %.2f, max %.2f nm' % (d[:,1].mean(), d[:,1].min(), d[:,1].max()))
print('   last 50%% mean %.2f nm' % d[t>=50,1].mean())
c=load('ncontacts.xvg')
print('contacts within cutoff: mean %.0f, min %.0f, max %.0f' % (c[:,1].mean(), c[:,1].min(), c[:,1].max()))
print('   last 50%% mean %.0f' % c[t>=50,1].mean())
p=load('c320_lig.xvg'); t2=p[:,0]
mind=p[:,1:].min(axis=1)
print('min dist C320 SG - ligand: mean %.2f nm, min %.2f, max %.2f' % (mind.mean(), mind.min(), mind.max()))
print('   last 50%% mean %.2f nm' % mind[t2>=50].mean())
s=load('sasa_prot.xvg')
print('protein SASA: mean %.1f nm^2, last-50%% %.1f' % (s[:,1].mean(), s[s[:,0]>=50,1].mean()))
