import numpy as np
import matplotlib.pyplot as plt
from util import utils

t08 = np.load('/home/pszypryt/sdss_data/20121208/Blue-Fit.npz')
t10 = np.load('/home/pszypryt/sdss_data/20121210/Blue10-Fit.npz')
t11 = np.load('/home/pszypryt/sdss_data/20121211/seq5Blue-Fit.npz')

params08 = t08['params']
params10 = t10['params']
params11 = t11['params']

jd08 = t08['jd']
jd10 = t10['jd']
jd11 = t11['jd']

params = np.vstack([params08,params10,params11])
jd = np.append(jd08,jd10)
jd = np.append(jd,jd11)

amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
jd2 = (jd/0.01966127)%1.
iPeriod = np.array(jd/0.01966127,dtype=np.int)
iPeriod -= iPeriod[0]

fig = plt.figure()
ax = fig.add_subplot(111)
curve = amps*widths**2
curve1 =np.append(curve[0:510], curve[780:])
curve =np.append(curve1[0:510]/np.average(curve1[0:510]), curve1[510:1170]/np.average(curve1[510:1170]))
curve = np.append(curve,curve1[1170:]/np.average(curve1[1170:]))
jd2=np.append(jd2[0:510],jd2[780:])


#curve /= np.median(curve)
#amps /= np.median(amps)

fwhm = 2*np.sqrt(2*np.log(2))*widths#pixels
fwhm = 0.5*fwhm #arcsec
fwhm = widths
#medFwhm = utils.median_filterNaN(fwhm,size=5)
#meanFwhm = utils.mean_filterNaN(fwhm,size=5)

meanXpos = utils.mean_filterNaN(xpos,size=7)
meanYpos = utils.mean_filterNaN(ypos,size=7)
#curve/=np.median(curve)
fwhm/=np.median(fwhm)
numbins=155
binwidth = 1/float(numbins)
average_array = []

for i in range(numbins):
    out_values = np.where(np.logical_and(jd2 >= i*binwidth,jd2 < (i+1)*binwidth))[0]
    iCurve = curve[out_values]
    iCurve = iCurve[iCurve != 0]
#    iCurve = iCurve[iCurve < 700]
    bin_average = np.median(iCurve)
    average_array.append(bin_average)
x=np.linspace(0,1,numbins)
ax.plot(x-0.54,average_array,'k.')
plt.xlabel('Phase')
plt.ylabel('Scaled Photon Count')

for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(20)

plt.show()

#np.savez('fix.npz',widths=medFwhm,x=meanXpos,y=meanYpos)
