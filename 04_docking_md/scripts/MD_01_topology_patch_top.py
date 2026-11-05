import re
p='topol.top'
s=open(p).read()
assert 'LIG_atomtypes' not in s
s=s.replace('#include "amber14sb.ff/forcefield.itp"',
            '#include "amber14sb.ff/forcefield.itp"\n\n; GAFF2 atom types for BDPrB (acpype)\n#include "LIG_atomtypes.itp"',1)
s=s.replace('; Include water topology','; Ligand topology\n#include "LIG.itp"\n\n; Include water topology',1)
if '[ molecules ]' in s:
    s=s.replace('Protein_chain_A     1','Protein_chain_A     1\nLIG                 1',1)
open(p,'w').write(s)
print('ok')
