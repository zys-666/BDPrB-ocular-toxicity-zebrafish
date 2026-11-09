p="redock_overlay_labeled3.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("set label_position, [0.0, 2.6, 0.0]","set label_position, [0.0, 2.4, 5.0]")
open(p,"w",encoding="utf-8").write(s)
print("patched")
