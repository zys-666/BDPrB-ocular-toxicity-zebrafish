# --- fig1/fig2: system colours ---
p="mmpbsa_figs.py"
s=open(p,encoding="utf-8").read()
s=s.split("# ---------------- 3. per-residue top 15")[0] + 'print("done")\n'
s=s.replace('NAVY, RED = "#1f3a73", "#7b1113"','NAVY, RED = "#82C0EA", "#F5C157"')
open(p,"w",encoding="utf-8").write(s)

# --- fig3: solid system colour per figure ---
p="mmpbsa_fig3.py"
s=open(p,encoding="utf-8").read()
s=s.replace('    cols = plt.cm.RdYlBu(np.linspace(0.05, 0.55, len(vals)))',
'''    col = "#82C0EA" if tag == "BDPrB" else "#F5C157"
    cols = [col] * len(vals)''')
open(p,"w",encoding="utf-8").write(s)
print("patched")
