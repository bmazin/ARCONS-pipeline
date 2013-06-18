import numpy as np
import matplotlib.pyplot as plt
from util import utils
from scipy import pi,sqrt,exp
from matplotlib.ticker import MultipleLocator


### This is a narrow-purpose code that was altered from the more general (and useful) plotFullData.py. It is used to compare a few different lightcurves from the same data. In this case, aperture vs psf fit (in addition to integration time 3 and 0.1).


BDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAllInt3newwvlcal.npz'
RDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAllInt.1newwvlcal.npz'
GDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8ImageStackAllInt3Aperture510.npz'
#RDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfRedInt3.npz'
#GDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfGreenInt3.npz'

intTime8 = 3
SintTime8 = .1


BDataFile8 = np.load(BDec8npzData)
Bparams8 = BDataFile8['params']
#jd8 = BDataFile8['jd']
jd8 = BMJD8=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3MJD_TDB.txt')
Bamps8 = Bparams8[:,1]
Bwidths8 = Bparams8[:,4]
Bcurve8 = 2*pi*Bamps8*Bwidths8**2/intTime8
hrs8= (jd8*24-int(jd8[0]*24))

RDataFile8 = np.load(RDec8npzData)
Rparams8 = RDataFile8['params']
Ramps8 = Rparams8[:,1]
Rwidths8 = Rparams8[:,4]
Rcurve8 = 2*pi*Ramps8*Rwidths8**2/SintTime8
#Sjd8 = RDataFile8['jd']
Sjd8 = BMJD8=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAllInt.1MJD_TDB.txt')
#Sjd8=MJD8 = np.array(Sjd8)-2400000.5
Shrs8= (Sjd8*24-int(Sjd8[0]*24))
#np.savez('/Scratch/dataProcessing/SDSS_J0926/Dec8AllInt.1lightcurve.npz',MJD_TDB = Sjd8,curve = Rcurve8)

GDataFile8 = np.load(GDec8npzData)
#Gparams8 = GDataFile8['params']
#Gamps8 = Gparams8[:,1]
#Gwidths8 = Gparams8[:,4]
#Gcurve8 = 2*pi*Gamps8*Gwidths8**2/intTime8
Gcurve8 = GDataFile8['lightcurve']


f = plt.figure()
plt.subplots_adjust(hspace=0.001)

ax1 = plt.subplot(111)
plt.title('Dec8 300-800 nm')
ax1.plot(Shrs8[:1300],Rcurve8[:1300],'r-',label= "psf fit 0.1 sec integration")
ax1.plot(Shrs8[1301:],Rcurve8[1301:],'r-')
ax1.plot(hrs8[:1300],Bcurve8[:1300],'b-',label= "psf fit 3 sec integration")
ax1.plot(hrs8[1301:],Bcurve8[1301:],'b-')
ax1.plot(hrs8[:1300],Gcurve8[:1300],'g-',label= "Aperture 3 sec integration, radius1 = 5 pix, radius2 = 10 pix")
ax1.plot(hrs8[1301:],Gcurve8[1301:],'g-')

plt.legend(loc="lower right")

plt.xlabel('Time(h) barycentered')

print np.average(Bwidths8),np.average(Bamps8)
print np.average(Rwidths8),np.average(Ramps8)

plt.show()
