import numpy as np
import matplotlib.pyplot as plt
from util import utils

<<<<<<< HEAD
t = np.load('/home/pszypryt/sdss_data/20121210/Blue10-Fit.npz')
=======
FileName = '/Scratch/dataProcessing/SDSS_J0926/20121208/ShortIntfitpsfBlue.npz'
NumFrames = 1700
IntTime = 3
FoldPeriod = 0.01966127 #This is in fractions of a day
t = np.load(FileName)
>>>>>>> c8b701cfcf2e94be0dbc57df049b7c2ec37a362f
params = t['params']
jd = t['jd']
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
jd2 = (jd/FoldPeriod)%1.

fig = plt.figure()
ax = fig.add_subplot(111)
#ax2 = fig.add_subplot(212)
curve = amps*widths**2
#curve /= np.median(curve)
#amps /= np.median(amps)
x = np.arange(0,NumFrames*IntTime,IntTime) # (0,930*10,10)
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
fwhm/=np.median(fwhm)
ax.plot(jd,curve,'k.')
<<<<<<< HEAD
=======
ax.set_title(FileName)
>>>>>>> c8b701cfcf2e94be0dbc57df049b7c2ec37a362f
#ax.plot(x,fwhm,'m')
#ax.plot(x,medFwhm,'k')
#ax.plot(x,meanFwhm,'b')
#plt.plot(x,amps)
#ax2.plot(jd[420:],xpos[420:],'.')
#ax2.plot(jd[420:],ypos[420:], '.')
#ax2.plot(jd,meanXpos)
#ax2.plot(jd,meanYpos)
plt.show()
#np.savez('fix.npz',widths=medFwhm,x=meanXpos,y=meanYpos)
