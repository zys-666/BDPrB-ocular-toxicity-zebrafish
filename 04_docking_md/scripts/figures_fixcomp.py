p="compose_labels.py"
s=open(p,encoding="utf-8").read()
s=s.replace('''for i in range(len(names)):
    for j in range(i + 1, len(names)):
        a, b = names[i], names[j]
        if abs(pts[a][0] - pts[b][0]) < 150 and abs(pts[a][1] - pts[b][1]) < 45:
            pts[a][1] -= 26
            pts[b][1] += 26''',
'''for _ in range(6):
    moved = False
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            if abs(pts[a][0] - pts[b][0]) < 175 and abs(pts[a][1] - pts[b][1]) < 52:
                cy = (pts[a][1] + pts[b][1]) / 2.0
                pts[a][1] = cy - 34
                pts[b][1] = cy + 34
                moved = True
    if not moved:
        break''')
open(p,"w",encoding="utf-8").write(s)
print("patched")
