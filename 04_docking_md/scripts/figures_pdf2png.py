import fitz, os
base=r"D:\user data\Desktop\轮胎\8.MD\mmpbsa"
for n in ["Graph1","Graph2"]:
    d=fitz.open(os.path.join(base,n+".pdf"))
    print(n, "pages:", d.page_count, "size:", d[0].rect)
    p=d[0].get_pixmap(dpi=150)
    out=os.path.join(r"D:\user data\Desktop\斑马鱼\5.dockMD\figures", n+"_ref.png")
    p.save(out)
    print("  saved", out, p.width, p.height)
