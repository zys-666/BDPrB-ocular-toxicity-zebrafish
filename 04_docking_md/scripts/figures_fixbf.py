p="hotspot_bfactor.py"
s=open(p,encoding="utf-8").read()
s=s.replace('f.write(l[:54] + "%6.2f" % b + l[60:].rstrip() + "\\n")',
            'f.write(l[:60] + "%6.2f" % b + l[66:].rstrip() + "\\n")')
open(p,"w",encoding="utf-8").write(s)
print("patched")
