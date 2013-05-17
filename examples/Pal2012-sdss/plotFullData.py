import numpy as np
import matplotlib.pyplot as plt
from util import utils
from scipy import pi,sqrt,exp
from matplotlib.ticker import MultipleLocator

BDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3.npz'
RDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfRedInt3.npz'
GDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfGreenInt3.npz'
R_Dec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfRed-Int3.npz'
I_Dec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfInfra-Int3.npz'
intTime8 = 3
#npzDataFile10 = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlue.npz'
BDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10.npz'
RDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfRedInt10.npz'
intTime10 = 10
#npzDataFile11 = '/Scratch/dataProcessing/SDSS_J0926/OlderData/Dec11fitpsfBlueUpdated.npz'
BDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10.npz'
RDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfRedInt10.npz'
intTime11 = 10

BDataFile8 = np.load(BDec8npzData)
Bparams8 = BDataFile8['params']
#jd8 = DataFile8['jd']
jd8 = BMJD8=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3MJD_TDB.txt')
#jd8 = np.array([line[1] for line in BMJD8])
#print jd8[0],jd8[len(jd8)-1]
#params[:,0] are the height offset of the gaussian (so the number of background pulses)
Bamps8 = Bparams8[:,1]
Bwidths8 = Bparams8[:,4]
#xpos8 = params8[:,2]
#ypos8 = params8[:,3]
Bcurve8 = 2*pi*Bamps8*Bwidths8**2/intTime8
hrs8= (jd8*24-int(jd8[0]*24))
#MJD8 = np.array(jd8)-2400000.5
#np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueMJDlong.txt',MJD8,fmt='%.30f')
#np.savez('/Scratch/dataProcessing/SDSS_J0926/Dec8BlueInt3lightcurve.npz',MJD_TDB = jd8,curve = curve8)

RDataFile8 = np.load(RDec8npzData)
Rparams8 = RDataFile8['params']
Ramps8 = Rparams8[:,1]
Rwidths8 = Rparams8[:,4]
Rcurve8 = 2*pi*Ramps8*Rwidths8**2/intTime8

GDataFile8 = np.load(GDec8npzData)
Gparams8 = GDataFile8['params']
Gamps8 = Gparams8[:,1]
Gwidths8 = Gparams8[:,4]
Gcurve8 = 2*pi*Gamps8*Gwidths8**2/intTime8

R_DataFile8 = np.load(R_Dec8npzData)
R_params8 = R_DataFile8['params']
R_amps8 = R_params8[:,1]
R_widths8 = R_params8[:,4]
R_curve8 = 2*pi*R_amps8*R_widths8**2/intTime8

I_DataFile8 = np.load(I_Dec8npzData)
I_params8 = I_DataFile8['params']
I_amps8 = I_params8[:,1]
I_widths8 = I_params8[:,4]
I_curve8 = 2*pi*I_amps8*I_widths8**2/intTime8

BDataFile10 = np.load(BDec10npzData)
Bparams10 = BDataFile10['params']
#jd10 = DataFile10['jd']
jd10 = BMJD10=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10MJD_TDB.txt')
#jd10 = np.array([line[1] for line in BMJD10])
Bamps10 = Bparams10[:,1]
Bwidths10 = Bparams10[:,4]
Bcurve10 = 2*pi*Bamps10*Bwidths10**2/intTime10
#MJD10 = np.array(jd10)-2400000.5
#np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10MJDlong.txt',MJD10,fmt='%.30f')

Bcurve10 = Bcurve10[421:]
jd10 = jd10[421:]
hrs10= (jd10*24-int(jd10[0]*24))
np.savez('/Scratch/dataProcessing/SDSS_J0926/Dec10BlueInt10lightcurve.npz',MJD_TDB = jd10,curve = Bcurve10)

RDataFile10 = np.load(RDec10npzData)
Rparams10 = RDataFile10['params']
Ramps10 = Rparams10[:,1]
Rwidths10 = Rparams10[:,4]
Rcurve10 = 2*pi*Ramps10*Rwidths10**2/intTime10
Rcurve10 = Rcurve10[421:]

BDataFile11 = np.load(BDec11npzData)
Bparams11 = BDataFile11['params']
#jd11 = DataFile11['jd']
jd11 = BMJD11=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10MJD_TDB.txt')
#jd11 = np.array([line[1] for line in BMJD11])
Bamps11 = Bparams11[:,1]
Bwidths11 = Bparams11[:,4]
Bcurve11 = 2*pi*Bamps11*Bwidths11**2/intTime11
hrs11= (jd11*24-int(jd11[0]*24))
#MJD11 = np.array(jd11)-2400000.5
#np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10MJDlong.txt',MJD11,fmt='%.30f')
#np.savez('/Scratch/dataProcessing/SDSS_J0926/Dec11BlueInt10lightcurve.npz',MJD_TDB = jd11,curve = curve11)

RDataFile11 = np.load(RDec11npzData)
Rparams11 = RDataFile11['params']
Ramps11 = Rparams11[:,1]
Rwidths11 = Rparams11[:,4]
Rcurve11 = 2*pi*Ramps11*Rwidths11**2/intTime11

f = plt.figure()
plt.subplots_adjust(hspace=0.001)


ax1 = plt.subplot(311)
ax1.plot(hrs8[:1300],Bcurve8[:1300],'b-',label="300 - 500 nm")
ax1.plot(hrs8[1301:],Bcurve8[1301:],'b-')
ax1.plot(hrs8[:1300],Rcurve8[:1300],'r-',label="500 - 700 nm")
ax1.plot(hrs8[1301:],Rcurve8[1301:],'r-')
ax1.plot(hrs8[:1300],Gcurve8[:1300],'g-',label="500 - 600 nm")
ax1.plot(hrs8[1301:],Gcurve8[1301:],'g-')
ax1.plot(hrs8[:1300],R_curve8[:1300],'Coral',label="600 - 700 nm")
ax1.plot(hrs8[1301:],R_curve8[1301:],'Coral')
ax1.plot(hrs8[:1300],I_curve8[:1300],'m-',label="700 - 900 nm")
ax1.plot(hrs8[1301:],I_curve8[1301:],'m-')
#plt.yticks(np.arange(700, 1700, 200))
plt.annotate('8th December 2012\nIntegration Time = 3 sec\nAverage Blue Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(Bwidths8),np.average(Bamps8)), xy=(0.7, 0.70), xycoords='axes fraction')
#plt.ylim(700,1600)
plt.ylim(100,1600)
plt.legend(loc="lower right")

ax2 = plt.subplot(312, sharex=ax1)
#ax2.plot(hrs10[:270],curve10[:270],'b-')
#ax2.plot(hrs10[271:420],curve10[271:420],'b-')
#ax2.plot(hrs10[421:],curve10[421:],'b-')
ax2.plot(hrs10,Bcurve10,'b-')
ax2.plot(hrs10,Rcurve10,'r-')
#plt.yticks(np.arange(0, 900, 200))
#plt.minorticks_on()
plt.annotate('10th December 2012\nIntegration Time = 10 sec\nAverage Blue Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(Bwidths10),np.average(Bamps10)), xy=(0.7, 0.70), xycoords='axes fraction')
#plt.ylim(0,900)
plt.ylim(0,500)
plt.ylabel('N photons Integrated Under Gaussian per sec')

ax3 = plt.subplot(313, sharex=ax1)
ax3.plot(hrs11,Bcurve11,'b-')
ax3.plot(hrs11,Rcurve11,'r-')
#plt.yticks(np.arange(0, 900, 200))
#plt.minorticks_on()
plt.annotate('11th December 2012\nIntegration Time = 10 sec\nAverage Blue Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(Bwidths11),np.average(Bamps11)), xy=(0.7, 0.70), xycoords='axes fraction')
#plt.ylim(0,900)
plt.ylim(0,500)
plt.xlabel('Time(h) barycentered')

xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
plt.setp(xticklabels, visible=False)

print np.average(Bwidths8),np.average(Bamps8)
print np.average(Bwidths10[421:]),np.average(Bamps10[421:])
print np.average(Bwidths11),np.average(Bamps11)

plt.show()
