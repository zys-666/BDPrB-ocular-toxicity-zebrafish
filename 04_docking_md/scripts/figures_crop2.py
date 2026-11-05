from PIL import Image
import os
p=r"D:\user data\Desktop\斑马鱼\13.图片\组图4.png"
im=Image.open(p); W,H=im.size
outdir=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\check"
crops={
 "cB":(int(W*0.30),0,int(W*0.55),int(H*0.13)),
 "cD2d":(int(W*0.10),int(H*0.10),int(W*0.32),int(H*0.24)),
 "cS":(int(W*0.74),int(H*0.63),W,int(H*0.80)),
 "cT":(0,int(H*0.80),W,H),
 "cMN":(int(W*0.55),int(H*0.40),W,int(H*0.56)),
}
for k,box in crops.items():
    c=im.crop(box)
    if c.width>1500: c=c.resize((1500,int(c.height*1500/c.width)), Image.LANCZOS)
    f=os.path.join(outdir,k+".png"); c.save(f); print(k,c.size)
