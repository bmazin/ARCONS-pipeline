import numpy as np
import matplotlib.pyplot as plt
from util import utils
from util.readDict import readDict
from scipy import pi

param = readDict()
param.read_from_file('0926params.dict')

FileName = param['npzLoadFitFile']
FramesPerFile = param['FramesPerFile']
print FramesPerFile
TotalNumFiles = param['TotalNumFiles']
NumFrames = FramesPerFile*TotalNumFiles
IntTime = param['integrationTime']

#FoldPeriod = 0.01966127 #This is in fractions of a day
t = np.load(FileName)
params = t['params']
jd = t['jd']
#params[:,0] are the height offset of the gaussian (so the number of background pulses)
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
#jd2 = (jd/FoldPeriod)%1.

fig = plt.figure()
ax = fig.add_subplot(111)
curve = 2*pi*amps*widths**2
#need to convert to Jy: curve*wvlrange*aveE/(1Hz*intTime*Area)

#fwhm = 2*np.sqrt(2*np.log(2))*widths#pixels
#fwhm = 0.5*fwhm #arcsec
#fwhm = widths

meanXpos = utils.mean_filterNaN(xpos,size=7)
meanYpos = utils.mean_filterNaN(ypos,size=7)
jd = jd[curve<=500]
curve = curve[curve<=500]
jd = jd[curve>=100]
curve = curve[curve>=100]

#curve/=np.median(curve)
#fwhm/=np.median(fwhm)
ax.plot(jd,curve,'k.')
ax.set_title(FileName)
plt.show()
