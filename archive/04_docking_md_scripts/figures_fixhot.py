p="hotspot_BDPrB.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("spectrum b, red_white_blue, hs and chain A, minimum=-8, maximum=2",
            "spectrum b, red_white_blue, hs and chain A, minimum=-8, maximum=8")
s=s.replace("set transparency, 0.25","set transparency, 0.35")
s=s.replace("zoom (hs and chain L), 16","zoom (hs and chain L), 26")
open(p,"w",encoding="utf-8").write(s)
print("patched")
