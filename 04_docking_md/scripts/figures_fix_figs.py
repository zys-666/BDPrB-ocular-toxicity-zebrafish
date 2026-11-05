p="md_figs.py"
s=open(p,encoding="utf-8").read()
s=s.replace('''for r in POCKET:
    ax.axvspan(r - 0.5, r + 0.5, color="#A5B59A", alpha=0.45, lw=0)''',
'''for lo, hi in [(136, 148), (186, 196), (260, 264), (312, 323), (473, 485)]:
    ax.axvspan(lo, hi, color="#A5B59A", alpha=0.35, lw=0)
ax.set_ylim(0, 6)
ax.text(505, 5.6, "C-terminal tail up to 17 A (cut)", fontsize=6, ha="right", color="0.4")''')
s=s.replace('ax[1].set_ylim(20, 32)','ax[1].set_ylim(0, 32)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
