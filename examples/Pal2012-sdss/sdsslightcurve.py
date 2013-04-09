import numpy as np
import matplotlib.pyplot as plt
from util import utils


FileName = '/Scratch/dataProcessing/SDSS_J0926/AllData/Dec11SIfitpsfRed.npz'
NumFrames = 2600
IntTime = 3

FoldPeriod = 0.01966127 #This is in fractions of a day
t = np.load(FileName)
params = t['params']
jd = t['jd']
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
jd2 = (jd/FoldPeriod)%1.

fig = plt.figure()
ax = fig.add_subplot(111)
curve = amps*widths**2
#curve /= np.median(curve)

fwhm = 2*np.sqrt(2*np.log(2))*widths#pixels
fwhm = 0.5*fwhm #arcsec
fwhm = widths

meanXpos = utils.mean_filterNaN(xpos,size=7)
meanYpos = utils.mean_filterNaN(ypos,size=7)

curve/=np.median(curve)
fwhm/=np.median(fwhm)
ax.plot(jd,curve,'k.')
ax.set_title(FileName)
plt.show()
