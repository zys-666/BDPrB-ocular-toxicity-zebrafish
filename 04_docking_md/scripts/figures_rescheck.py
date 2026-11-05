seen={}
for l in open("../dock/docking_results/complex_BDPrB.pdb"):
    if l.startswith("ATOM") and l[21]=="A":
        r=int(l[22:26]); seen.setdefault(r, l[17:20].strip())
for r in [138,142,146,187,188,191,192,195,262,286,287,314,319,320,321,342,475,476,477,478,479,483]:
    print(r, seen.get(r,"-"))
