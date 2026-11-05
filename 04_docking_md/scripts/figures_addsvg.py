p="fig56.py"
s=open(p,encoding="utf-8").read()
if "fig5_contacts.svg" not in s:
    s=s.replace('plt.savefig("fig5_contacts.png", dpi=300)',
                'plt.savefig("fig5_contacts.png", dpi=300)\nplt.savefig("fig5_contacts.svg")')
open(p,"w",encoding="utf-8").write(s)
print("patched")
