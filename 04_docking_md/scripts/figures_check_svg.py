import re, glob
for f in ["fig1_rmsd_rg.svg","fig2_ligand.svg","fig3_rmsf.svg","fig4_pocket_bars.svg","fig5_contacts.svg","complex_BDPrB.svg","complex_KyA33.svg"]:
    c=open(f,encoding="utf-8",errors="replace").read()
    nt=len(re.findall(r"<text", c)); npth=len(re.findall(r"<path", c))
    tnr = "Times New Roman" in c or "TimesNewRoman" in c
    grid = "grid" in c.lower()
    print("%-22s text=%3d path=%3d  TimesNR=%s  'grid'=%s" % (f, nt, npth, tnr, grid))
