import numpy as np
for tag,fn in [("BDPrB","s9_heat_BDPrB.npz"),("KyA33","s9_heat_KyA33.npz")]:
    d=np.load(fn); m=d["mat"]; r=d["res"]
    per_frame=m.sum(0)
    print(tag, "res", m.shape, "per-frame contacts min/med/max:",
          int(per_frame.min()), int(np.median(per_frame)), int(per_frame.max()))
    top=np.argsort(per_frame)[-5:]
    print("   frames with most contacts:", [(int(i), int(per_frame[i])) for i in top])
    print("   residues ever contacting:", int((m.sum(1)>0).sum()))
