from PIL import Image, ImageDraw, ImageFont
import os
base=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures"
a=Image.open(os.path.join(base,"bead_overview_KyA33.png")).convert("RGB")
b=Image.open(os.path.join(base,"redock_overlay_toon.png")).convert("RGB")
H=900
a=a.resize((int(a.width*H/a.height),H), Image.LANCZOS)
b=b.resize((int(b.width*H/b.height),H), Image.LANCZOS)
pad=30; label=90
W=a.width+b.width+pad*3
canvas=Image.new("RGB",(W,H+pad*2),(255,255,255))
canvas.paste(a,(pad,pad))
canvas.paste(b,(pad*2+a.width,pad))
d=ImageDraw.Draw(canvas)
try:
    f=ImageFont.truetype("arialbd.ttf", 56)
except:
    f=ImageFont.load_default()
d.text((pad+10,pad+10), "A", fill=(40,40,40), font=f)
d.text((pad*2+a.width+10,pad+10), "B", fill=(40,40,40), font=f)
out=os.path.join(base,"fig_redock_2panel.png")
canvas.save(out, dpi=(300,300))
print("saved", out, canvas.size)
