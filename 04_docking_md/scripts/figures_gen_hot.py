# -*- coding: utf-8 -*-
"""Generate toon-style hot-spot pmls (base + one-label probes)."""
import os

RES = [(138, "VAL-138"), (188, "PHE-188"), (195, "TRP-195"),
       (314, "PHE-314"), (477, "MET-477")]

HEAD = """load ../dock/docking_results/complex_BDPrB.pdb, cplx
hide everything
show cartoon, cplx and chain A
color 0xCFEAF6, cplx and chain A
set cartoon_fancy_helices, 1
set cartoon_flat_sheets, 1

show sticks, cplx and chain A and resi 138+188+195+314+477
color 0xA9A2D6, cplx and chain A and resi 138+188+195+314+477
set stick_radius, 0.16, cplx and chain A and resi 138+188+195+314+477
hide (hydro and cplx and chain A)

show sticks, cplx and chain L
color 0xC0392B, cplx and chain L
set stick_radius, 0.30, cplx and chain L
hide (hydro and cplx and chain L)

set ray_trace_mode, 1
set ray_trace_color, grey25
set ray_trace_gain, 0.08
set ray_shadows, 0
set specular, 0
set reflect, 0
set ambient, 0.62
set direct, 0.42
set light_count, 1
set depth_cue, 0
set antialias, 2
bg_color white
set ray_opaque_background, 1

orient cplx and chain A
zoom (cplx and chain L), 7.5
viewport 3200, 2560
"""

open("hot_base.pml", "w", encoding="utf-8").write(
    HEAD + "png hot_base.png, dpi=600, ray=1\nquit\n")

PROBE_HEAD = HEAD.replace("set ray_trace_mode, 1", "set ray_trace_mode, 0")
PROBE_HEAD = PROBE_HEAD.replace("color 0xCFEAF6, cplx and chain A", "color white, cplx and chain A")
PROBE_HEAD = PROBE_HEAD.replace("color 0xC0392B, cplx and chain L", "color white, cplx and chain L")
PROBE_HEAD = PROBE_HEAD.replace("color 0xA9A2D6, cplx and chain A and resi 138+188+195+314+477",
                                "color white, cplx and chain A and resi 138+188+195+314+477")

for r, t in RES:
    s = PROBE_HEAD
    s += "set label_size, 30\nset label_color, black\nset label_position, [0.0, 2.0, 0.0]\n"
    s += 'label cplx and chain A and resi %d and name CA, "%s"\n' % (r, t)
    s += "png hot_probe_%d.png, width=3200, height=2560, ray=0\nquit\n" % r
    open("hot_probe_%d.pml" % r, "w", encoding="utf-8").write(s)
print("generated", len(RES) + 1, "pml files")
