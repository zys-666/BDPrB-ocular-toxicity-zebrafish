from PIL import Image, ImageDraw, ImageFont
a = Image.open("bead_overview_KyA33.png").convert("RGB")
b = Image.open("redock_overlay_toon.png").convert("RGB")
H = 900
a = a.resize((int(a.width * H / a.height), H), Image.LANCZOS)
b = b.resize((int(b.width * H / b.height), H), Image.LANCZOS)
pad = 30
W = a.width + b.width + pad * 3
canvas = Image.new("RGB", (W, H + pad * 2), (255, 255, 255))
canvas.paste(a, (pad, pad))
canvas.paste(b, (pad * 2 + a.width, pad))
d = ImageDraw.Draw(canvas)
try:
    f = ImageFont.truetype("arialbd.ttf", 56)
except Exception:
    f = ImageFont.load_default()
d.text((pad + 10, pad + 10), "A", fill=(40, 40, 40), font=f)
d.text((pad * 2 + a.width + 10, pad + 10), "B", fill=(40, 40, 40), font=f)
canvas.save("fig_redock_2panel.png", dpi=(300, 300))
print("saved", canvas.size)
