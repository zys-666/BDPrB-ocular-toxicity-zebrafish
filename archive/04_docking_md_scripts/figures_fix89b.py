p="s89_plot.py"
s=open(p,encoding="utf-8").read()
# crop whitespace around the structure image
s=s.replace('img = Image.open(os.path.join(BASE, "5.dockMD", "figures", "s8_structure.png")).convert("RGB")',
'''img = Image.open(os.path.join(BASE, "5.dockMD", "figures", "s8_structure.png")).convert("RGB")
import numpy as _np
_a = _np.array(img.convert("L"))
_ys, _xs = _np.where(_a < 250)
_pad = 25
img = img.crop((max(0, _xs.min() - _pad), max(0, _ys.min() - _pad),
                min(img.width, _xs.max() + _pad), min(img.height, _ys.max() + _pad)))''')
# more room for colorbars
s=s.replace('plt.subplots_adjust(left=0.09, right=0.97, top=0.93, bottom=0.07)',
            'plt.subplots_adjust(left=0.09, right=0.93, top=0.93, bottom=0.07)')
open(p,"w",encoding="utf-8").write(s)
print("patched")
