import numpy as np, glob
for f in ["chk.pdb","chk2.pdb"]:
    frames=[]; cur=None
    for l in open(f):
        if l.startswith("MODEL"): cur={"p":[],"l":[]}
        elif l.startswith("ENDMDL"): frames.append(cur); cur=None
        elif l.startswith(("ATOM","HETATM")) and cur is not None:
            x=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
            (cur["l"] if l[17:20].strip()=="KYA" else cur["p"]).append(x)
    ds=[np.linalg.norm(np.array(fr["l"])[:,None]-np.array(fr["p"])[None],axis=-1).min() for fr in frames]
    print(f, len(frames), "frames, min distances:", [round(d,2) for d in ds])
