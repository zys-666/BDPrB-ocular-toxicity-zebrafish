p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\build_overlay_complex.py"
s=open(p,encoding='utf-8').read()
s=s.replace('''                l[:6].strip().ljust(6), n, name, " ", l[17:20].strip(), chain,''',
'''                ("HETATM" if chain in ("C","D") else "ATOM"), n, name, " ",
                ("LIG" if chain in ("C","D") else l[17:20].strip()), chain,''')
open(p,'w',encoding='utf-8').write(s)
print("patched")
