import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy import pi


def Ephem(x,Period,Pdot):
#    return 53795.9455191+Period*(x)+.5*(Pdot)*(x)**2

#    Period = ReportedPeriod
    jd0 = 53795.9455191
    return jd0+Period*x/(1.-Pdot*x)

def linearEphem(x,Period):
#    return 53795.9455191+Period*(x)+.5*(Pdot)*(x)**2

#    Period = ReportedPeriod
    jd0 = 53795.9455191
    return jd0+Period*x

def fitEphemeris(x,data,sigma):
    p0=(ReportedPeriod,-10**-7)
    popt, pcov = optimize.curve_fit(Ephem, x, data, p0=p0,sigma=sigma)
    print p0
    return popt

def fitLinearEphemeris(x,data,sigma):
    p0=(ReportedPeriod)
    popt, pcov = optimize.curve_fit(linearEphem, x, data, p0=p0,sigma=sigma)
    print p0
    return popt

ReportedPeriod = 0.01966127289

DataFileBlue_ = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataBlue-fitresults.npz')
paramsBlue_ = DataFileBlue_['FittedParams']
EclipseTimesBlue_ = DataFileBlue_['EclipseTimes']
DataFileBlue = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataBluefitresults.npz')
paramsBlue = DataFileBlue['FittedParams']
EclipseTimesBlue = DataFileBlue['EclipseTimes']
DataFileRed = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataRedfitresults.npz')
paramsRed = DataFileRed['FittedParams']
EclipseTimesRed = DataFileRed['EclipseTimes']
DataFileMult = np.load('/Scratch/dataProcessing/SDSS_J0926/Dec8AddColorsfitresults.npz')
paramsMult = DataFileMult['FittedParams']
EclipseTimesMult = DataFileMult['EclipseTimes']
EclipseNumMult = [0,1,2,0,1,2,0,1,2]

OldData = np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/CopperwheatResults.dat')
OldEclipseNum = np.array([line[0] for line in OldData])
OldEclipseTime = np.array([line[1] for line in OldData])
OldEclipseError = np.array([line[2] for line in OldData])

gap1=100
gap2=47
Eclipse=[0,1,2,gap1+2,gap1+3,gap1+4,gap1+4+gap2,gap1+5+gap2,gap1+6+gap2,gap1+7+gap2,gap1+8+gap2]

EclipseTimes = [EclipseTimesBlue,EclipseTimesRed,EclipseTimesMult,EclipseTimesBlue_]
EclipseNum = [Eclipse,Eclipse,EclipseNumMult,Eclipse]
FixedEclipseTimes = [[],[],[],[]]
FixedEclipseNum = [[],[],[],[]]

for ind in range(len(EclipseTimes)):
    for eclipse in range(len(EclipseTimes[ind])):
        if np.isnan(EclipseTimes[ind][eclipse])==False:
            FixedEclipseTimes[ind].append(EclipseTimes[ind][eclipse])
            FixedEclipseNum[ind].append(EclipseNum[ind][eclipse])

x0=125860
FixedEclipseNum = np.array([np.array(line)+x0 for line in FixedEclipseNum])

TotalEclipseNum = np.hstack((OldEclipseNum,FixedEclipseNum[3]))
TotalEclipseTime = np.hstack((OldEclipseTime,FixedEclipseTimes[3]))

pOpt = fitEphemeris(TotalEclipseNum,TotalEclipseTime,None)
print pOpt

linearPeriod = fitLinearEphemeris(TotalEclipseNum,TotalEclipseTime,None)

Std = np.std(FixedEclipseTimes[3]-Ephem(FixedEclipseNum[3],ReportedPeriod,0))*24*60*60
ave = np.average(FixedEclipseTimes[3]-Ephem(FixedEclipseNum[3],ReportedPeriod,0))*24*60*60

fig=plt.figure()
plt.subplots_adjust(hspace=.5,wspace=.6)
ax=fig.add_subplot(321)
ax.plot(FixedEclipseNum[3],(FixedEclipseTimes[3]-Ephem(FixedEclipseNum[3],ReportedPeriod,0))*24*60*60,'b*',label="Blue' 400-550 nm with new wvl sol")
ax.errorbar(np.average(FixedEclipseNum[3]),ave,yerr=Std,fmt='bo',label="Blue' average = %f with std = %f"%(ave,Std))
plt.legend(loc='upper left',prop={'size':'small'})
ax2=fig.add_subplot(322)
ax2.plot(FixedEclipseNum[0],(FixedEclipseTimes[0]-Ephem(FixedEclipseNum[0],ReportedPeriod,0))*24*60*60,'bo',label="Blue 300-500 nm")
ax2.plot(FixedEclipseNum[1],(FixedEclipseTimes[1]-Ephem(FixedEclipseNum[1],ReportedPeriod,0))*24*60*60,'r^',label="Red 500-700 nm")
plt.legend(loc='upper left',prop={'size':'small'})
ax3=fig.add_subplot(323)
ax3.plot(FixedEclipseNum[0],(FixedEclipseTimes[0]-Ephem(FixedEclipseNum[0],ReportedPeriod,0))*24*60*60,'bo',label="Blue 300-500 nm")
ax3.plot(FixedEclipseNum[2][:3],(FixedEclipseTimes[2][:3]-Ephem(FixedEclipseNum[2][:3],ReportedPeriod,0))*24*60*60,'g>',label="Green 300-400 nm")
ax3.plot(FixedEclipseNum[2][3:],(FixedEclipseTimes[2][3:]-Ephem(FixedEclipseNum[2][3:],ReportedPeriod,0))*24*60*60,'m<',label="Red' 400-500 nm")
plt.legend(loc='upper left',prop={'size':'small'})

ax4=fig.add_subplot(324)
ax4.set_title('Divergence from Reported Ephemeris: Data - Equation',size='small')

ax4.errorbar(OldEclipseNum,(OldEclipseTime-Ephem(OldEclipseNum,ReportedPeriod,0))*24*60*60, yerr=OldEclipseError, fmt='g^',label="Copperwheat Data")
ax4.plot(FixedEclipseNum[3],(FixedEclipseTimes[3]-Ephem(FixedEclipseNum[3],ReportedPeriod,0))*24*60*60,'b*',label="Our Data from 400-550nm")
ax4.set_ylabel('O-C (s)')
ax4.set_xlabel('Cycle')
plt.legend(loc='upper left',prop={'size':'small'})
ax5=fig.add_subplot(325)
ax5.errorbar(OldEclipseNum,(OldEclipseTime-Ephem(OldEclipseNum,pOpt[0],pOpt[1]))*24*60*60, yerr=OldEclipseError, fmt='g^',label="Copperwheat Data")
ax5.plot(FixedEclipseNum[3],(FixedEclipseTimes[3]-Ephem(FixedEclipseNum[3],pOpt[0],pOpt[1]))*24*60*60,'b*',label="Our Data from 400-550nm")
ax5.set_xlabel('Cycle')
ax5.set_ylabel('O-C (s)')
ax5.set_title(' Reported Period = 0.01966127289, Fitted Period = %.10f, Fitted Pdot = %.15f\nPeriod Difference (Reported-Fitted) = %.10f days or %.10f seconds\nFit does not take error into account\n\nDivergence from Fitted Ephemeris'%(pOpt[0],pOpt[1],ReportedPeriod-pOpt[0],(ReportedPeriod-pOpt[0])*24*60*60),size='small')
plt.legend(loc='upper left',prop={'size':'small'})

ax6=fig.add_subplot(326)
ax6.set_title(' Reported Period = 0.01966127289, Fitted Period = %.10f\nPeriod Difference (Reported-Fitted) = %.10f days or %.10f seconds\nFit does not take error into account\n\nDivergence from Linear Fitted Ephemeris: Data - Equation'%(linearPeriod,ReportedPeriod-linearPeriod,(ReportedPeriod-linearPeriod)*24*60*60),size='small')

ax6.errorbar(OldEclipseNum,(OldEclipseTime-Ephem(OldEclipseNum,linearPeriod,0))*24*60*60, yerr=OldEclipseError, fmt='g^',label="Copperwheat Data")
ax6.plot(FixedEclipseNum[3],(FixedEclipseTimes[3]-Ephem(FixedEclipseNum[3],linearPeriod,0))*24*60*60,'b*',label="Our Data from 400-550nm")
ax4.set_ylabel('O-C (s)')
ax4.set_xlabel('Cycle')
plt.legend(loc='upper left',prop={'size':'small'})

plt.show()




