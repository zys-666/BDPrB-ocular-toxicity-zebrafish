p=r"D:\user data\Desktop\斑马鱼\5.dockMD\figures\make_illustrate_pdb.py"
s=open(p,encoding="utf-8").read()
s=s.replace("ks = sorted(set(zc) & set(hc))\nP = np.array([hc[k] for k in ks]); Q = np.array([zc[k] for k in ks])",
            "OFF = 23  # 9R6N PDB numbering = zebrafish/UniProt numbering - 23\nks = sorted(k for k in hc if (k + OFF) in zc)\nP = np.array([hc[k] for k in ks]); Q = np.array([zc[k + OFF] for k in ks])")
open(p,'w',encoding="utf-8").write(s)
print("patched")
