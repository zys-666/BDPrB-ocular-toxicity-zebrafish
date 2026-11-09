p="s89_plot.py"
s=open(p,encoding="utf-8").read()
s=s.replace('ax.text(13, 88, "1-25\\nremoved", ha="center", fontsize=8, color="#8E2B22", fontweight="bold")\n', "")
s=s.replace('ax.text(13, 15, "mean 41.0", ha="center", fontsize=7.5, color="#8E2B22")',
'''ax.annotate("N-terminal 1-25 removed\\n(mean pLDDT 41.0)", xy=(22, 42), xytext=(80, 22),
            fontsize=8, color="#8E2B22",
            arrowprops=dict(arrowstyle="->", color="#8E2B22", lw=0.9))''')
s=s.replace('''ax2.text(0.02, 0.03, "N-terminal 1-25 (red) reaches the catalytic cavity;\\nCys320 (orange sphere)",
         transform=ax2.transAxes, fontsize=7.5, color="0.2", va="bottom")''',
'''ax2.text(0.5, -0.02, "N-terminal 1-25 (red) reaches the catalytic cavity; Cys320 (orange sphere)",
         transform=ax2.transAxes, fontsize=7.5, color="0.2", ha="center", va="top")''')
s=s.replace('plt.subplots_adjust(left=0.08, right=0.99, top=0.90, bottom=0.13)',
            'plt.subplots_adjust(left=0.08, right=0.99, top=0.90, bottom=0.16)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
