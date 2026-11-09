p="s89_plot.py"
s=open(p,encoding="utf-8").read()
s=s.replace('plt.subplots_adjust(left=0.09, right=0.93, top=0.93, bottom=0.07)',
            'plt.subplots_adjust(left=0.09, right=0.90, top=0.93, bottom=0.07)')
s=s.replace('cb.set_label("contact", fontsize=7.5)','cb.set_label("contact", fontsize=7.5, labelpad=1)')
open(p,"w",encoding="utf-8").write(s)
