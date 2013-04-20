import numpy as np
import matplotlib.pyplot as plt
from util import utils

cutoff = 75
TemplateFileName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTemplateUpdated100.npz'
TruncTempFileName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTruncTemplateUpdated100.npz'

#t08 = np.load('/home/pszypryt/sdss_data/20121208/Blue-Fit.npz')
#t10 = np.load('/home/pszypryt/sdss_data/20121210/Blue10-Fit.npz')
#t11 = np.load('/home/pszypryt/sdss_data/20121211/seq5Blue-Fit.npz')
t10 = np.load('/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8fitpsfBlueUpdated.npz')

#params08 = t08['params']
params10 = t10['params']
#params11 = t11['params']

#jd08 = t08['jd']
jd10 = t10['jd']
#jd11 = t11['jd']

#params = np.vstack([params08,params10,params11])
#jd = np.append(jd08,jd10)
#jd = np.append(jd,jd11)

params = params10
jd = jd10

period = 0.01966127

amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
jd2 = (jd/period)%1.

#find number of points per period, so bin density is the same as data density (1 data point = 1 bin)
#maxcheck = np.max(jd2[:600])
#print maxcheck
#nbins = np.argmax(jd2[:600])
#print nbins
#maxcheck1 = np.max(jd2[601:1200])
#print maxcheck1
#nbins1 = np.argmax(jd2[601:1200])
#nbins1 = nbins1+600-nbins
#print nbins1
#maxcheck2 = np.max(jd2[1201:])
#print maxcheck2
#nbins2 = np.argmax(jd2[1201:])
#print nbins2+1200-nbins-nbins1

#plt.plot(jd2)
#plt.show()

iPeriod = np.array(jd/period,dtype=np.int)
iPeriod -= iPeriod[0]

fig = plt.figure()
ax = fig.add_subplot(211)
curve1 = amps*widths**2
curve = np.append(curve1[0:411]/np.average(curve1[0:411]), curve1[411:971]/np.average(curve1[411:971])) 
curve = np.append(curve, curve1[971:]/np.average(curve1[971:]))
#curve1 =np.append(curve[0:510], curve[780:])
#curve =np.append(curve1[0:510]/np.average(curve1[0:510]), curve1[510:1170]/np.average(curve1[510:1170]))
#curve = np.append(curve,curve1[1170:]/np.average(curve1[1170:]))
#jd2=np.append(jd2[0:510],jd2[780:])


#curve /= np.median(curve)
#amps /= np.median(amps)

#fwhm = 2*np.sqrt(2*np.log(2))*widths#pixels
#fwhm = 0.5*fwhm #arcsec
fwhm = widths
#medFwhm = utils.median_filterNaN(fwhm,size=5)
#meanFwhm = utils.mean_filterNaN(fwhm,size=5)

#meanXpos = utils.mean_filterNaN(xpos,size=7)
#meanYpos = utils.mean_filterNaN(ypos,size=7)
#curve/=np.median(curve)
#fwhm/=np.median(fwhm)
#numbins=560 one bin per point
numbins = 230
binwidth = 1/float(numbins)
jdbin = period/float(numbins)
average_array = []

for i in range(numbins):
    out_values = np.where(np.logical_and(jd2 >= i*binwidth,jd2 < (i+1)*binwidth))[0]
    iCurve = curve[out_values]
    iCurve = iCurve[iCurve != 0]
#    iCurve = iCurve[iCurve < 700]
    bin_average = np.mean(iCurve)
    average_array.append(bin_average)

#jd=np.arange(0,1,binwidth)*period
jd3 = np.arange(0,period,jdbin)
print len(jd), len(jd3)
print jd[105], jd3[105]



np.savez(TemplateFileName,template=average_array,jd=jd3) 
MinIndex= np.argmin(average_array)

#ax.plot(x-0.54,average_array,'k.')
ax.plot(jd3,average_array,'k.')
ax.set_xlabel('Phase')
ax.set_ylabel('Scaled Photon Count')
plt.title(TruncTempFileName)
ax2=fig.add_subplot(212)
ax2.plot(jd3[MinIndex-cutoff:MinIndex+cutoff],average_array[MinIndex-cutoff:MinIndex+cutoff], 'k.')
ax2.set_xlabel('Phase')
ax2.set_ylabel('Scaled Photon Count')

for item in ([ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels() + [ax2.xaxis.label, ax2.yaxis.label] +
             ax2.get_xticklabels() + ax2.get_yticklabels()):
    item.set_fontsize(20)

np.savez(TruncTempFileName,template=average_array[MinIndex-cutoff:MinIndex+cutoff],jd=jd3[MinIndex-cutoff:MinIndex+cutoff]) 

plt.show()


#np.savez('fix.npz',widths=medFwhm,x=meanXpos,y=meanYpos)
