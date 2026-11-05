p=r"D:\user data\Desktop\斑马鱼\5.dockMD\controls\compare_md.py"
s=open(p,encoding='utf-8').read()
s=s.replace('"BDPrB": os.path.join(root, "MD", "03_production"),',
            '"BDPrB": os.path.join(root, "MD", "04_analysis"),')
open(p,'w',encoding='utf-8').write(s)
print('ok')
