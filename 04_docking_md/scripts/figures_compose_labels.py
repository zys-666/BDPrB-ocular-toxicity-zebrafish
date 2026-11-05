# -*- coding: utf-8 -*-
"""Composite residue labels on top of the clean render (labels never occluded)."""
import numpy as np
from PIL import Image, ImageDraw, ImageFont

items = [(188, "PHE188"), (192, "MET192"), (195, "TRP195"), (314, "PHE314"),
         (319, "CYS319"), (320, "CYS320"), (321, "THR321"), (475, "ASN475")]

pos = {}
for r, t in items:
    im = Image.open("probe_%d.png" % r).convert("L")
    a = np.array(im)
    ys, xs = np.where(a < 128)
    if len(xs) == 0:
        print("no text found for", t)
        continue
    pos[t] = (int(xs.mean()), int(ys.mean()), int(xs.min()), int(xs.max()), int(ys.min()), int(ys.max()))
    print(t, pos[t])

base = Image.open("base_clean.png").convert("RGB")
d = ImageDraw.Draw(base)
try:
    font = ImageFont.truetype(r"C:\Windows\Fonts\timesbd.ttf", 40)
except Exception:
    font = ImageFont.load_default()

# nudge labels that would collide
pts = {t: [p[0], p[1]] for t, p in pos.items()}
names = list(pts)
for _ in range(6):
    moved = False
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            if abs(pts[a][0] - pts[b][0]) < 350 and abs(pts[a][1] - pts[b][1]) < 104:
                cy = (pts[a][1] + pts[b][1]) / 2.0
                pts[a][1] = cy - 68
                pts[b][1] = cy + 68
                moved = True
    if not moved:
        break

for t in names:
    x, y = pts[t]
    bb = d.textbbox((0, 0), t, font=font)
    w, h = bb[2] - bb[0], bb[3] - bb[1]
    d.text((x - w / 2, y - h / 2 - bb[1]), t, font=font, fill=(25, 25, 25),
           stroke_width=6, stroke_fill=(255, 255, 255))

base.save("redock_overlay_labeled.png", dpi=(600, 600))
print("saved redock_overlay_labeled.png", base.size)
