from PIL import Image
p=r"D:\user data\Desktop\斑马鱼\13.图片\组图4.png"
im=Image.open(p); print(im.size, im.mode)
W,H=im.size
crops={
 "c1_top":(0,0,W,int(H*0.22)),
 "c2_row2":(0,int(H*0.20),int(W*0.62),int(H*0.40)),
 "c3_2d1":(int(W*0.28),int(H*0.20),int(W*0.62),int(H*0.40)),
 "c4_2d2":(int(W*0.60),int(H*0.06),W,int(H*0.40)),
}
import os
outdir=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\check"
os.makedirs(outdir,exist_ok=True)
for k,box in crops.items():
    c=im.crop(box)
    if c.width>1600: c=c.resize((1600,int(c.height*1600/c.width)), Image.LANCZOS)
    f=os.path.join(outdir,k+".png"); c.save(f); print(k, c.size, f)
