head = open("probe.pml",encoding="utf-8").read().split("python")[0]
items = [(188,"PHE188"),(192,"MET192"),(195,"TRP195"),(314,"PHE314"),
         (319,"CYS319"),(320,"CYS320"),(321,"THR321"),(475,"ASN475")]
for r,t in items:
    s = head + 'label cplx and chain A and resi %d and name CA, "%s"\n' % (r,t)
    s += 'png probe_%d.png, width=1800, height=1400, ray=0\nquit\n' % r
    open("probe_%d.pml"%r,"w",encoding="utf-8").write(s)
print("wrote 8 probe pml")
