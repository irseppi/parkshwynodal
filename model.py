import numpy as np
import matplotlib.pyplot as plt

fnot = [37*2, 56*2, 73*2, 110*2, 146*2, 165*2, 182*2, 218*2, 238*2, 256*2, 275*2]

tprime0 = 107
tpr = np.arange(0, 241, 1)
c = 343
v0 = 100
l = 2700
plt.figure(figsize=(10, 6))
for f0 in fnot:
	ft = []
	for tprime in tpr:
		ft0p = f0*1/(1+(v0/c)*(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))/(np.sqrt(l**2+(v0*((tprime - tprime0)- np.sqrt((tprime-tprime0)**2-(1-v0**2/c**2)*((tprime-tprime0)**2-l**2/c**2)))/(1-v0**2/c**2))**2)))
			
		ft.append(ft0p)
	plt.plot(tpr, ft, 'g', linewidth=0.5)

plt.axvline(x=tprime0, c = 'g', ls = '--', label='Estimated t0: '+str(tprime0)+' sec')

plt.ylabel('Frequency (Hz)')
plt.title("Forward Model: t'= "+str(tprime0)+' sec, v0 = '+str(int(v0))+' m/s, l = '+str(l)+' m, \n' + 'f0 = '+str(fnot)+' Hz', fontsize='small')

plt.xlim(0, 240)

plt.show()
					
