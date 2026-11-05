p="hotspot_BDPrB.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("spectrum b, red_white_blue, hs and chain A, minimum=-8, maximum=8",
            "spectrum b, red_white, hs and chain A, minimum=-7, maximum=0")
s=s.replace("set transparency, 0.35","set transparency, 0.10")
open(p,"w",encoding="utf-8").write(s)
print("patched")
