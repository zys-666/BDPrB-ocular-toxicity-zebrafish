import glob, re
# --- redock overlay: 600 dpi = 3600x2800 px, protein #CFEAF6
for p in ["base_clean.pml"] + glob.glob("probe_*.pml"):
    s=open(p,encoding="utf-8").read()
    s=s.replace("color 0xFBD9C9","color 0xCFEAF6")
    s=s.replace("viewport 1800, 1400","viewport 3600, 2800")
    s=s.replace("png base_clean.png, dpi=300","png base_clean.png, dpi=600")
    s=s.replace("dpi=300, ray=1","dpi=600, ray=1")
    if "png probe_" in s:
        s=re.sub(r"width=1800, height=1400","width=3600, height=2800",s)
    open(p,"w",encoding="utf-8").write(s)
# label size scales with pixels
s=open("compose_labels.py",encoding="utf-8").read()
s=s.replace('"timesbd.ttf", 40','"timesbd.ttf", 80').replace("dpi=(300, 300)","dpi=(600, 600)")
s=s.replace("stroke_width=3","stroke_width=6")
s=s.replace("< 175","< 350").replace("< 52","< 104").replace("- 34","- 68").replace("+ 34","+ 68")
open("compose_labels.py","w",encoding="utf-8").write(s)
print("patched redock")
