import numpy as np
import matplotlib.pyplot as plt
from util import utils

#This program folds, bins, and averages the data to create a template. It can be fit using fitTemplate.py and can be used in a convolution using ConvolutionWithTemp.py. Make sure that the numbins used in this program is the same value used in ConvolutionWithTemp so the data and template have the same point/time ratio.

cutoff = 75
TemplateFileName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTemplateUpdated100.npz'
TruncTempFileName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTruncTemplateUpdated100.npz'

#t08 = np.load('/home/pszypryt/sdss_data/20121208/Blue-Fit.npz')
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

iPeriod = np.array(jd/period,dtype=np.int)
iPeriod -= iPeriod[0]

fig = plt.figure()
ax = fig.add_subplot(211)
curve = amps*widths**2
curve /= np.average(curve)

fwhm = widths

numbins=560 #one bin per point
#numbins = 230
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

jd3 = np.arange(0,period,jdbin)

np.savez(TemplateFileName,template=average_array,jd=jd3) 
MinIndex= np.argmin(average_array)

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
