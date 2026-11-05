import re
for p in ["md_figs.py","fig56.py"]:
    s=open(p,encoding="utf-8").read()
    # rcParams: Times New Roman + text (not paths)
    head = '''import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
'''
    repl = '''import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = ["Times New Roman", "DejaVu Serif"]
plt.rcParams["mathtext.fontset"] = "stix"
'''
    if "svg.fonttype" not in s:
        s = s.replace(head, repl, 1)
    # kill grids
    s = re.sub(r"\.grid\(alpha=[^)]*\)", ".grid(False)", s)
    s = s.replace(".grid(alpha=0.25, lw=0.5, axis=\"y\")", ".grid(False)")
    s = s.replace(".grid(alpha=0.25, lw=0.5)", ".grid(False)")
    # superscript in fig4
    s = s.replace('"Pocket SASA (nm$^2$)"', '"Pocket SASA (nm²)"')
    open(p,"w",encoding="utf-8").write(s)
    print("patched", p, "svg.fonttype" in s, ".grid(False)" in s)
