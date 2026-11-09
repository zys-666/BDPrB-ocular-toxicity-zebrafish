p="s89_plot.py"
s=open(p,encoding="utf-8").read()
s=s.replace('"N-terminal 1-25 (red) reaches the catalytic cavity; Cys320 (orange sphere)"',
            '"N-terminal 1-25 (red) reaches the catalytic cavity;\\nCys320 shown as orange sphere"')
open(p,"w",encoding="utf-8").write(s)
print("patched")
