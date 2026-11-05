p='extract.py'
s=open(p,encoding='utf-8').read()
s=s.replace("""    for a in prot:
        f.write('ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s\\n' % (
            0, a['atom'], a['comp'], int(a['seq']), a['x'], a['y'], a['z'], a['b'],""",
"""    for i, a in enumerate(prot, 1):
        f.write('ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00%6.2f          %2s\\n' % (
            i, a['atom'], a['comp'], int(a['seq']), a['x'], a['y'], a['z'], a['b'],""")
open(p,'w',encoding='utf-8').write(s)
print('patched')
