p="s89_plot.py"
s=open(p,encoding="utf-8").read()
s=s.replace('    sm = np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 1, mat)',
            '    sm = np.apply_along_axis(lambda m: np.convolve(m, k, mode="same"), 1, mat).T')
s=s.replace('axA.set_ylim(0, 100)','axA.set_ylim(0, 5)')
s=s.replace('axA.text(0.5, 0.62, "BDPrB: 0%\\n(no donor, ether O only)", transform=axA.transAxes,',
            'axA.text(0.5, 0.72, "BDPrB: 0%\\n(no donor; ether O only)", transform=axA.transAxes,')
open(p,"w",encoding="utf-8").write(s)
print("patched")
