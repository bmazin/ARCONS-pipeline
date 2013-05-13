import numpy as np
import matplotlib.pyplot as plt
from util import utils
from scipy import pi,sqrt,exp
from matplotlib.ticker import MultipleLocator

BDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3.npz'
RDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfRedInt3.npz'
intTime8 = 3
#npzDataFile10 = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlue.npz'
BDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10.npz'
RDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfRedInt10.npz'
intTime10 = 10
#npzDataFile11 = '/Scratch/dataProcessing/SDSS_J0926/OlderData/Dec11fitpsfBlueUpdated.npz'
BDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10.npz'
intTime11 = 10

DataFile8 = np.load(BDec8npzData)
params8 = DataFile8['params']
jd8 = DataFile8['jd']
#BMJD8=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueBMJD.txt')
#jd8 = np.array([line[1] for line in BMJD8])
#print jd8[0],jd8[len(jd8)-1]
#params[:,0] are the height offset of the gaussian (so the number of background pulses)
amps8 = params8[:,1]
widths8 = params8[:,4]
#xpos8 = params8[:,2]
#ypos8 = params8[:,3]
curve8 = 2*pi*amps8*widths8**2/intTime8
hrs8= (jd8*24-int(jd8[0]*24))
MJD8 = np.array(jd8)-2400000.5
np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueMJDlong.txt',MJD8,fmt='%.30f')

DataFile10 = np.load(BDec10npzData)
params10 = DataFile10['params']
jd10 = DataFile10['jd']
#BMJD10=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10BMJD.txt')
#jd10 = np.array([line[1] for line in BMJD10])
amps10 = params10[:,1]
widths10 = params10[:,4]
curve10 = 2*pi*amps10*widths10**2/intTime10
MJD10 = np.array(jd10)-2400000.5
np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10MJDlong.txt',MJD10,fmt='%.30f')

curve10 = curve10[421:]
jd10 = jd10[421:]
hrs10= (jd10*24-int(jd10[0]*24))

DataFile11 = np.load(BDec11npzData)
params11 = DataFile11['params']
jd11 = DataFile11['jd']
#BMJD11=np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10BMJD.txt')
#jd11 = np.array([line[1] for line in BMJD11])
amps11 = params11[:,1]
widths11 = params11[:,4]
curve11 = 2*pi*amps11*widths11**2/intTime11
hrs11= (jd11*24-int(jd11[0]*24))
MJD11 = np.array(jd11)-2400000.5
np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10MJDlong.txt',MJD11,fmt='%.30f')

f = plt.figure()
plt.subplots_adjust(hspace=0.001)


ax1 = plt.subplot(311)
ax1.plot(hrs8[:1300],curve8[:1300],'b-')
ax1.plot(hrs8[1301:],curve8[1301:],'b-')
plt.yticks(np.arange(700, 1700, 200))
#plt.minorticks_on()
#minorLocator   = MultipleLocator(1)
#ax1.xaxis.set_minor_locator(minorLocator)
plt.annotate('8th December 2012\nIntegration Time = 3 sec\nAverage Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(widths8),np.average(amps8)), xy=(0.7, 0.70), xycoords='axes fraction')
plt.ylim(700,1600)

ax2 = plt.subplot(312, sharex=ax1)
#ax2.plot(hrs10[:270],curve10[:270],'b-')
#ax2.plot(hrs10[271:420],curve10[271:420],'b-')
#ax2.plot(hrs10[421:],curve10[421:],'b-')
ax2.plot(hrs10,curve10,'b-')
plt.yticks(np.arange(0, 900, 200))
#plt.minorticks_on()
plt.annotate('10th December 2012\nIntegration Time = 10 sec\nAverage Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(widths10),np.average(amps10)), xy=(0.7, 0.70), xycoords='axes fraction')
plt.ylim(0,900)
plt.ylabel('N photons Integrated Under Gaussian per sec')

ax3 = plt.subplot(313, sharex=ax1)
ax3.plot(hrs11,curve11)
plt.yticks(np.arange(0, 900, 200))
#plt.minorticks_on()
plt.annotate('11th December 2012\nIntegration Time = 10 sec\nAverage Gaussian Parameters:\nWidth = %.2f, Height = %.2f'%(np.average(widths11),np.average(amps11)), xy=(0.7, 0.70), xycoords='axes fraction')
plt.ylim(0,900)
plt.xlabel('Time(h) barycentered')

xticklabels = ax1.get_xticklabels()+ax2.get_xticklabels()
plt.setp(xticklabels, visible=False)

print np.average(widths8),np.average(amps8)
print np.average(widths10[421:]),np.average(amps10[421:])
print np.average(widths11),np.average(amps11)

plt.show()
