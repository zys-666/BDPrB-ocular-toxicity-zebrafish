p="md_figs.py"
s=open(p,encoding="utf-8").read()
s=s.replace('ax[1].plot(g[:, 0], g[:, 1] * 10, lw=0.7, color=C[k], label=k)',
            'ax[1].plot(g[:, 0] / 1000.0, g[:, 1] * 10, lw=0.7, color=C[k], label=k)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
