import re
# --- hotspot ---
p="gen_hot.py"
s=open(p,encoding="utf-8").read()
s=s.replace("color 0xFBD9C9","color 0xCFEAF6")
s=s.replace("viewport 1600, 1280","viewport 3200, 2560")
s=s.replace("png hot_base.png, dpi=300","png hot_base.png, dpi=600")
s=s.replace("width=1600, height=1280","width=3200, height=2560")
open(p,"w",encoding="utf-8").write(s)

p="compose_hot3.py"
s=open(p,encoding="utf-8").read()
s=s.replace('"timesbd.ttf", 40','"timesbd.ttf", 80')
s=s.replace("stroke_width=3","stroke_width=6")
s=s.replace("<210","<420").replace("<52","<104").replace("-34","-68").replace("+34","+68")
s=s.replace("dpi=(300,300)","dpi=(600,600)")
open(p,"w",encoding="utf-8").write(s)

# --- snapshots ---
p="snap5.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("0xFBD9C9","0xCFEAF6")
s=s.replace("viewport 800, 800","viewport 1600, 1600")
s=s.replace("dpi=200","dpi=600")
open(p,"w",encoding="utf-8").write(s)
print("patched hotspot + snapshots")
