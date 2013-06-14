import numpy as np
import matplotlib.pyplot as plt
from util import utils
import mpfit
import scipy.optimize as optimize

#This program slides a template across the data to try to improve the data quality. This convultion was not very successful. The program bins the data if neccessary so the template and data have the same point/time ratio. "Trunc" means the template does not cover the entire period, just the area around the eclipse. The number in the template filenames indicate the number of bins per period. dataTemplate.py makes a template by binning and averaging the data. fitTemplate.py makes a template by fitting a skewed gaussian to the "data template" that results from templateMaker.py

FileName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8fitpsfBlueUpdated.npz'
NumFrames = 1700
IntTime = 3
TimeCutOff = 1300

template = np.load('/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTruncTemplateFitted560.npz')
template3 = np.load('/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8BlueTruncTemplateUpdated560.npz')

tempValues = template['template']
tempjd = template['jd']

tempValues2 = template2['template']
tempjd2 = template2['jd']

tempValues3 = template3['template']
tempjd3 = template3['jd']

t = np.load(FileName)

params = t['params']
jd = t['jd']

period = 0.01966127

amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]

NumFrames = TimeCutOff
jd = jd[:TimeCutOff]
amps = amps[:TimeCutOff]
widths = widths[:TimeCutOff]
xpos = xpos[:TimeCutOff]
ypos = ypos[:TimeCutOff]

fig = plt.figure()
ax = fig.add_subplot(311)
curve1 = amps*widths**2
curve = curve/(np.average(curve))
#curve = np.append(curve1[0:411]/np.average(curve1[0:411]), curve1[411:971]/np.average(curve1[411:971])) 
#curve = np.append(curve, curve1[971:]/np.average(curve1[971:]))

numbins = 155 #per eclipse period!
jdbin = period/float(numbins)
Totalnumbins = int((jd[NumFrames-1]-jd[0])/jdbin)

#to bin the DATA
average_array = []

for i in range(Totalnumbins):
    out_values = np.where(np.logical_and(jd >= i*jdbin+jd[0],jd < (i+1)*jdbin+jd[0]))[0]
    iCurve = curve[out_values]
    iCurve = iCurve[iCurve != 0]
#    iCurve = iCurve[iCurve < 700]
    bin_average = np.mean(iCurve)
    average_array.append(bin_average)

#jd=np.arange(0,1,binwidth)*period
jd3 = np.arange(jd[0]+jdbin/2.,jd[NumFrames-1]-jdbin/2.,jdbin)

print len(average_array),len(jd3),Totalnumbins

ax.plot(jd,curve,'g')
#ax.plot(jd3,average_array,'r')
ax.set_title(FileName)
#plt.show()

ax2 = fig.add_subplot(312)

MinIndex= np.argmin(tempValues)
tempValues4 = tempValues3[MinIndex-20:MinIndex+20]
tempjd4 = tempjd3[MinIndex-20:MinIndex+20]

filtered = np.correlate(tempValues,curve,mode='same')[::-1]
filtered2 = np.correlate(tempValues2,curve,mode='same')[::-1]
filtered3 = np.correlate(tempValues3,curve,mode='same')[::-1]
filtered4 = np.correlate(tempValues4,curve,mode='same')[::-1]

ax.plot(jd,filtered/np.average(filtered),'r')
ax.plot(jd,filtered4/np.average(filtered4),'k')
ax.plot(jd,filtered3/np.average(filtered3),'b')
#ax2.plot(tempjd2,tempValues2,'.r')
ax2.plot(jd,filtered/np.average(filtered),'r')
ax2.plot(jd,filtered4/np.average(filtered4)+.1,'k')
ax2.plot(jd,filtered3/np.average(filtered3)+.2,'b')

ax3 = fig.add_subplot(313)
ax3.plot(tempjd,tempValues,'.r',label="skew-gauss fitted")
ax3.plot(tempjd4,tempValues4,'ko', label="short not fitted")
ax3.plot(tempjd3,tempValues3,'.b', label="not fitted")
#ax.plot(jd[0]+tempjd2,tempValues2,'.k')
#ax.plot(jd[0]+.01+tempjd,tempValues,'.r')
plt.legend(loc=4)

plt.show()



