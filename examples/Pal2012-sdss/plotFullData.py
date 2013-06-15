import numpy as np
import matplotlib.pyplot as plt
from util import utils
from scipy import pi,sqrt,exp
from matplotlib.ticker import MultipleLocator
import matplotlib

### This program is used to create a plot showing the light curve of multiple different light bands on different days. It also has the ability to save a txt file of the data times, which are then sent to tempo2 to be barycentered and then ConvertTCBtoTDB.py to switch it over to TDB. It can also save an npz file with the lightcurve only.

#Dec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAllInt3newwvlcal.npz'
BDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3.npz'
#RDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAllInt.1newwvlcal.npz'
#GDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8ImageStackAllInt3Aperture510.npz'

intTime8 = 3
#SintTime8 = .1
#npzDataFile10 = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlue.npz'
BDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfAllInt10newwvlcal.npz'
#BDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10.npz'
#RDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfRedInt10newwvlcal.npz'
intTime10 = 10
#BDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/OlderData/Dec11fitpsfBlueUpdated.npz'
BDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10.npz'
#BDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfAllInt10newwvlcal.npz'
#RDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfRedInt10.npz'
intTime11 = 10

BDataFile8 = np.load(BDec8npzData)
Bparams8 = BDataFile8['params']
#jd8 = BDataFile8['jd']
jd8 = BMJD8=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3MJD_TDB.txt')
#params[:,0] are the height offset of the gaussian (so the number of background pulses)
Bamps8 = Bparams8[:,1]
Bwidths8 = Bparams8[:,4]
Bcurve8 = 2*pi*Bamps8*Bwidths8**2/intTime8
Bcurve8 /= np.average(Bcurve8)
hrs8= (jd8*24-int(jd8[0]*24))
#jd8=MJD8 = np.array(jd8)-2400000.5
#np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAll.1MJDlong.txt',MJD8,fmt='%.30f')
#np.savez('/Scratch/dataProcessing/SDSS_J0926/Dec8AllInt.1lightcurve.npz',MJD_TDB = jd8,curve = Bcurve8)

#RDataFile8 = np.load(RDec8npzData)
#Rparams8 = RDataFile8['params']
#Ramps8 = Rparams8[:,1]
#Rwidths8 = Rparams8[:,4]
#Rcurve8 = 2*pi*Ramps8*Rwidths8**2/intTime8
#Sjd8 = RDataFile8['jd']
#Sjd8 = BMJD8=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAllInt.1MJD_TDB.txt')
#Sjd8=MJD8 = np.array(Sjd8)-2400000.5
#Shrs8= (Sjd8*24-int(Sjd8[0]*24))
#np.savez('/Scratch/dataProcessing/SDSS_J0926/Dec8AllInt.1lightcurve.npz',MJD_TDB = Sjd8,curve = Rcurve8)

#GDataFile8 = np.load(GDec8npzData)
#Gparams8 = GDataFile8['params']
#Gamps8 = Gparams8[:,1]
#Gwidths8 = Gparams8[:,4]
#Gcurve8 = 2*pi*Gamps8*Gwidths8**2/intTime8
#Gcurve8 = GDataFile8['lightcurve']

BDataFile10 = np.load(BDec10npzData)
Bparams10 = BDataFile10['params']
#jd10 = DataFile10['jd']
jd10 = BMJD10=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10MJD_TDB.txt')
Bamps10 = Bparams10[:,1]
Bwidths10 = Bparams10[:,4]
Bcurve10 = 2*pi*Bamps10*Bwidths10**2/intTime10
#Bcurve10 /= np.average(Bcurve10[421:])
#MJD10 = np.array(jd10)-2400000.5
#np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10MJDlong.txt',MJD10,fmt='%.30f')

Bcurve10 = Bcurve10[421:]
Bcurve10 /= np.average(Bcurve10)
jd10 = jd10[421:]
hrs10= (jd10*24-int(jd10[0]*24))
#np.savez('/Scratch/dataProcessing/SDSS_J0926/Dec10BlueInt10lightcurve.npz',MJD_TDB = jd10,curve = Bcurve10)

#RDataFile10 = np.load(RDec10npzData)
#Rparams10 = RDataFile10['params']
#Ramps10 = Rparams10[:,1]
#Rwidths10 = Rparams10[:,4]
#Rcurve10 = 2*pi*Ramps10*Rwidths10**2/intTime10
#Rcurve10 = Rcurve10[421:]

BDataFile11 = np.load(BDec11npzData)
Bparams11 = BDataFile11['params']
#jd11 = DataFile11['jd']
jd11 = BMJD11=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10MJD_TDB.txt')
Bamps11 = Bparams11[:,1]
Bwidths11 = Bparams11[:,4]
Bcurve11 = 2*pi*Bamps11*Bwidths11**2/intTime11
Bcurve11 /= np.average(Bcurve11)
hrs11= (jd11*24-int(jd11[0]*24))
#MJD11 = np.array(jd11)-2400000.5
#np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10MJDlong.txt',MJD11,fmt='%.30f')
#np.savez('/Scratch/dataProcessing/SDSS_J0926/Dec11BlueInt10lightcurve.npz',MJD_TDB = jd11,curve = curve11)

#RDataFile11 = np.load(RDec11npzData)
#Rparams11 = RDataFile11['params']
#Ramps11 = Rparams11[:,1]
#Rwidths11 = Rparams11[:,4]
#Rcurve11 = 2*pi*Ramps11*Rwidths11**2/intTime11



fig = plt.figure()
font = {'size'   : 20}
matplotlib.rc('font', **font)
plt.subplots_adjust(hspace=0.001)

ax1 = plt.subplot(211)
plt.ylim(0.5,1.5)
plt.yticks([0.5,1.0])
ax1.plot(hrs8[:1300],Bcurve8[:1300],'b.')
ax1.plot(hrs8[1301:],Bcurve8[1301:],'b.')

#ax1.plot(hrs8[:1300],Gcurve8[:1300],'g-',label= "Aperture 3 sec integration, radius1 = 5 pix, radius2 = 10 pix") #"500 - 600 nm")
#ax1.plot(hrs8[1301:],Gcurve8[1301:],'g-')
#ax1.plot(hrs8[:1300],Rcurve8[:1300],'Coral',label="600 - 700 nm")
#ax1.plot(hrs8[1301:],Rcurve8[1301:],'Coral')

#plt.yticks(np.arange(700, 1700, 200))
plt.annotate('8th December 2012\nIntegration Time = 3 sec', xy=(0.05, 0.75), xycoords='axes fraction')#\nAverage Blue Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(Bwidths8),np.average(Bamps8)),)
#plt.ylim(100,1600)
#plt.legend(loc="lower right")

#ax2 = plt.subplot(312, sharex=ax1, sharey=ax1)

#ax2.plot(hrs10[:270],Bcurve10[:270],'b-')
#ax2.plot(hrs10[271:420],Bcurve10[271:420],'b-')
#ax2.plot(hrs10[421:],Bcurve10[421:],'b-')

#plt.yticks(np.arange(0, 900, 200))
#plt.minorticks_on()
#plt.annotate('10th December 2012\nIntegration Time = 10 sec', xy=(0.7, 0.10), xycoords='axes fraction')#\nAverage Blue Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(Bwidths10),np.average(Bamps10)))
#plt.ylim(0,1.5)
#plt.yticks([0,0.5,1.0])


ax3 = plt.subplot(212, sharex=ax1, sharey=ax1)
ax3.xaxis.labelpad = 20
ax3.yaxis.labelpad = 20
ax3.plot(hrs11,Bcurve11,'b.')
#ax3.plot(hrs11,Rcurve11,'r-')
#plt.yticks([0,0.5,1.0])
#plt.minorticks_on()
plt.annotate('11th December 2012\nIntegration Time = 10 sec', xy=(0.05, 0.75), xycoords='axes fraction')#\nAverage Blue Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(Bwidths11),np.average(Bamps11)))
plt.minorticks_on()
ax3.set_xlabel('Time (hours)')
ax3.set_ylabel('                           Normalized Flux')
xticklabels = ax1.get_xticklabels()
plt.setp(xticklabels, visible=False)

#print np.average(Bwidths8),np.average(Bamps8)
#print np.average(Bwidths10[421:]),np.average(Bamps10[421:])
#print np.average(Bwidths11),np.average(Bamps11)

plt.ylim(0.5,1.5)
plt.yticks([0.5,1.0])

plt.show()
