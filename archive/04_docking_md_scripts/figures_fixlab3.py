p="redock_overlay_labeled3.pml"
s=open(p,encoding="utf-8").read()
s=s.replace("set label_position, [0.0, 2.4, 20.0]","set label_position, [0.0, 2.0, 0.0]")
s=s.replace("set ray_trace_mode, 1","set ray_trace_mode, 0")
open(p,"w",encoding="utf-8").write(s)
print("patched")
