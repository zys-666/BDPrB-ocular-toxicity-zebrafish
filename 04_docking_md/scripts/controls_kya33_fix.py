p=r'D:\user data\Desktop\斑马鱼\5.dockMD\controls\kya33\map_pose.py'
s=open(p,encoding='utf-8').read()
s=s.replace("""def coords(path):
    r = []
    for l in open(path):
        if l.startswith(('ATOM', 'HETATM')):
            r.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(r)""",
"""def coords(path):
    r = []
    for l in open(path):
        if l.startswith(('ATOM', 'HETATM')):
            el = l[77:79].strip()
            if el in ('H', 'HD'):
                continue
            r.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(r)""")
s=s.replace("pose = pose[:n_at]","pose = pose[:n_at]")
open(p,'w',encoding='utf-8').write(s)
print('patched')
