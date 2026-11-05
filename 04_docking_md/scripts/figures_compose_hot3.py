import numpy as np
from PIL import Image, ImageDraw, ImageFont

items = [(138,"VAL-138"),(188,"PHE-188"),(195,"TRP-195"),(314,"PHE-314"),(477,"MET-477")]
blank = np.array(Image.open("hot_blank.png").convert("L")).astype(int)
pts = {}
for r,t in items:
    a = np.array(Image.open("hot_probe_%d.png"%r).convert("L")).astype(int)
    diff = (blank - a) > 40
    ys,xs = np.where(diff)
    pts[t] = [float(xs.mean()), float(ys.mean())]
    print(t, [round(v) for v in pts[t]], "pixels", len(xs))

base = Image.open("hot_base.png").convert("RGB")
d = ImageDraw.Draw(base)
try:
    font = ImageFont.truetype(r"C:\Windows\Fonts\timesbd.ttf", 40)
except Exception:
    font = ImageFont.load_default()
names = list(pts)
for _ in range(8):
    moved=False
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            a,b = names[i],names[j]
            if abs(pts[a][0]-pts[b][0])<420 and abs(pts[a][1]-pts[b][1])<104:
                cy=(pts[a][1]+pts[b][1])/2.0
                pts[a][1]=cy-68; pts[b][1]=cy+68; moved=True
    if not moved: break
for t in names:
    x,y = pts[t]
    bb = d.textbbox((0,0),t,font=font)
    w,h = bb[2]-bb[0], bb[3]-bb[1]
    d.text((x-w/2, y-h/2-bb[1]), t, font=font, fill=(20,20,20), stroke_width=6, stroke_fill=(255,255,255))
base.save("hotspot_toon.png", dpi=(600,600))
print("saved", base.size)
