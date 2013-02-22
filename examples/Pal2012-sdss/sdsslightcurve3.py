import numpy as np
import matplotlib.pyplot as plt
from util import utils

t = np.load('/home/pszypryt/sdss_data/20121211/seq5Blue-Fit.npz')
params = t['params']
jd = t['jd']
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
jd2 = (jd/0.01966127)%1.
iPeriod = np.array(jd/0.01966127,dtype=np.int)
iPeriod -= iPeriod[0]
jd_0 = jd2[iPeriod == 0]
jd_1 = jd2[iPeriod == 1]
jd_2 = jd2[iPeriod == 2]
jd_3 = jd2[iPeriod == 3]

fig = plt.figure()
ax = fig.add_subplot(111)
curve = amps*widths**2
#curve /= np.median(curve)
#amps /= np.median(amps)
x = np.arange(0,930*10,10)
#xpos/=np.median(xpos)
#ypos/=np.median(ypos)

fwhm = 2*np.sqrt(2*np.log(2))*widths#pixels
fwhm = 0.5*fwhm #arcsec
fwhm = widths
#medFwhm = utils.median_filterNaN(fwhm,size=5)
#meanFwhm = utils.mean_filterNaN(fwhm,size=5)

meanXpos = utils.mean_filterNaN(xpos,size=7)
meanYpos = utils.mean_filterNaN(ypos,size=7)

curve/=np.median(curve)
curve_0 = curve[iPeriod == 0]
curve_1 = curve[iPeriod == 1]
curve_2 = curve[iPeriod == 2]
curve_3 = curve[iPeriod == 3]
fwhm/=np.median(fwhm)
ax.plot(jd_0,curve_0,'c.')
ax.plot(jd_1,curve_1,'y.')
ax.plot(jd_2,curve_2,'b.')
ax.plot(jd_3,curve_3,'r.')
plt.show()
#np.savez('fix.npz',widths=medFwhm,x=meanXpos,y=meanYpos)
