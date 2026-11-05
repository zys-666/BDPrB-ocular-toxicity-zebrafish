import re
lines=open('LIG.itp').read().splitlines()
i=lines.index('[ atoms ]')
itp=[]
for l in lines[i+1:]:
    if l.strip().startswith('['): break
    p=l.split()
    if len(p)>=5 and not l.strip().startswith(';'):
        itp.append((int(p[0]), p[4], float(p[6])))
print('itp atoms:', len(itp), 'net ligand charge: %.4f' % sum(a[2] for a in itp))

g=open('solv_ions.gro').read().splitlines()
gr=[]
for l in g[2+7571:2+7571+50]:
    gr.append((int(l[15:20]), l[10:15].strip(), l[5:10].strip()))
ok = all(itp[k][0]==gr[k][0] and itp[k][1]==gr[k][2] for k in range(50))
print('itp/gro ligand atom order consistent:', ok)
for k in range(5):
    print('  ', itp[k], gr[k])
if not ok:
    for k in range(50):
        if itp[k][1]!=gr[k][2]: print('  MISMATCH', k, itp[k], gr[k])
