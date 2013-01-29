import numpy as np
import matplotlib.pyplot as plt
from util import utils

t = np.load('nlttFit5.npz')
t2 = np.load('nlttFit.npz')
params = t['params']
params2 = t2['params']
errors = t['errors']
jd = t['jd']

transitJD = 2456272.62693
chisq = t['chisqs'][:,0]
dof = t['chisqs'][:,1]
reducedChisq = chisq/dof
amps = params[:,1]
amps2 = params2[:,1]
widths = params[:,4]
widthErrs = errors[:,4]
widths2 = params2[:,4]
xpos = params[:,2]
xposErrs = errors[:,2]
ypos = params[:,3]
yposErrs = errors[:,3]
xpos2 = params2[:,2]
ypos2 = params2[:,3]

fig = plt.figure()
ax = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
curve = amps*widths**2
curve2 = amps2*widths2**2
#curve /= np.median(curve)
#amps /= np.median(amps)
#x = np.arange(0,720*15,15)
x = np.arange(0,360*30,30)
x2 = np.arange(0,360*30,30)
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
curve2/=np.median(curve2)
print np.shape(curve),np.shape(x)
ax.plot(x,curve,'k')
ax.plot(x[range(65,71)],curve[range(65,71)],'r')
ax.plot(x,widths,'r')
ax.plot(x,widths+widthErrs,'r')
ax.plot(x,widths-widthErrs,'r')
ax.plot(x2,widths2,'m')

#ax.plot(x,medFwhm,'k')
#ax.plot(x,meanFwhm,'b')
#plt.plot(x,amps)
ax2.plot(x,xpos)
ax2.plot(x,xpos+xposErrs)
ax2.plot(x,xpos-xposErrs)
ax2.plot(x2,xpos2)
ax2.plot(x,ypos)
ax2.plot(x,ypos+yposErrs)
ax2.plot(x,ypos-yposErrs)
ax2.plot(x2,ypos2)
#ax2.plot(x,meanXpos)
#ax2.plot(x,meanYpos)
ax3.plot(x,jd)
plt.show()
#np.savez('fix.npz',widths=medFwhm,x=meanXpos,y=meanYpos)
