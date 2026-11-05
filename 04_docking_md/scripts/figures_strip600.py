from PIL import Image, ImageDraw, ImageFont
files=[("snap/pocket_0.png","0 ns"),("snap/pocket_25.png","25 ns"),("snap/pocket_50.png","50 ns"),
       ("snap/pocket_75.png","75 ns"),("snap/pocket_100.png","100 ns")]
imgs=[Image.open(f).convert("RGB") for f,_ in files]
w,h=imgs[0].size
pad,lab=20,92
canvas=Image.new("RGB",(w*5+pad*6,h+lab+pad*2),(255,255,255))
d=ImageDraw.Draw(canvas)
try: fnt=ImageFont.truetype("arialbd.ttf",68)
except Exception: fnt=ImageFont.load_default()
for i,(im,(_,t)) in enumerate(zip(imgs,files)):
    x=pad+i*(w+pad)
    canvas.paste(im,(x,lab))
    d.text((x+24,20),t,fill=(40,40,40),font=fnt)
canvas.save("fig6_snapshots.png",dpi=(600,600))
print("saved fig6",canvas.size)
