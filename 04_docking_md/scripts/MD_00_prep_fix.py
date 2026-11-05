p='make_complex.py'
s=open(p,encoding='utf-8').read()
s=s.replace("out.append('%4d %4d 0 0 0' % (22, len(bonds)))","nb_written = sum(1 for b in bonds if int(b[1])<=22 and int(b[2])<=22)\nout.append('%4d %4d 0 0 0' % (22, nb_written))\nprint('heavy bonds:', nb_written)")
open(p,'w',encoding='utf-8').write(s)
