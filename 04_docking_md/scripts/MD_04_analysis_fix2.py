p='validate_structure.py'
s=open(p,encoding='utf-8').read()
s=s.replace("res.append((l[5:10].strip(), int(l[15:20]), l[10:15].strip(),",
            "res.append((l[5:10].strip(), int(l[0:5]), l[10:15].strip(),")
open(p,'w',encoding='utf-8').write(s)
print('fixed')
