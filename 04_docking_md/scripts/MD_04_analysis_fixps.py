p='run_analysis.ps1'
s=open(p,encoding='utf-8').read()
s=s.replace('@("Temperature","Pressure","Density","Potential") + "" | Set-Content s8.txt -Encoding ASCII',
            '@("Temperature","Pressure","Density","Potential","") | Set-Content s8.txt -Encoding ASCII')
open(p,'w',encoding='utf-8').write(s)
print('fixed')
