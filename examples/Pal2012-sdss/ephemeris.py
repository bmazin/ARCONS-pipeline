import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy import pi
import matplotlib

### np.polyfit is better at optimizing than scipy.optimize.curvefit
### there is a optimizing bug/confusion. Using PdotLimit.py which finds the optimum Pdot by finding the minimum reduduced chi-squared I found fDot=-1.2623365827*10**-9,Pdot=4.88*10**-13 is the best with chi-sq= 1.09. Using curvefit, fitted fDot=-3.20872712971*10**-10,Pdot=1.24*10**-13, chi-sq=1.15. Using polyfit, fDot=-1.25949536823*10**-9, Pdot=4.87*10**-13, chi-sq=1.18.
### so I don't know which is best. I'm going to have the program stick with the polyfit for now... 

#This program loads the eclipse times saved from DwarfFit.py and the data from Copperwheat and calculates and plots the observed - calculated eclipse times for each eclipse using copperwheat's ephemeris, fitted linear ephemeris, fitted quadratic ephemeris as the models for the calculated eclipse time. The fits are weighted by the error of the measurements which is the std of each set of data. The reduced chi-squared is given for each model.

#Load up the right fit results, check that Eclipse/EclipseShort (line 68/69) has the eclipse numbers in the right form. Alter EclipseTimes, EclipseNum to match what you want, and make sure you change lines 78-94 if you are using more than one fits file. The rest should be fine unchanged.

def Ephem(t,jd0,f,fDot):
    return f*(t-jd0) + 0.5 * fDot * (t-jd0)**2

def FixedEphem(t,jd0,f):
#    jd0 = 53795.9455191
    fDot=-1.2623365827*10**-9
    return f*(t-jd0) + 0.5 * fDot * (t-jd0)**2

def fitModel(t):
    return iNumber - Ephem(t,jd0,f,fDot)

def fitFixedEphem(t,data,sigma):
    p0=(ReportedOffset,1/ReportedPeriod)
    popt, pcov = optimize.curve_fit(FixedEphem, t, data, p0=p0,sigma=sigma)
    return popt

def fitEphem(t,data,sigma):
    p0=(ReportedOffset,1/ReportedPeriod,-10**-13/ReportedPeriod**2)
    popt, pcov = optimize.curve_fit(Ephem, t, data, p0=p0,sigma=sigma)
    return popt

def RedChiSq(diff,sigma,DoF):
    return np.sum([diff[ind]**2/(1.*sigma[ind]**2*DoF) for ind in range(len(diff))]) 
#    DoF = N(num of points)-v(num variables)-1

def linearEphem(x,jd0,Period):
    return jd0+Period*x

ReportedPeriod = 0.01966127289
ReportedOffset = 53795.9455191

DataFileAll = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataAllwvlfitresults.npz')
paramsAll = DataFileAll['FittedParams']
EclipseTimesAll = DataFileAll['EclipseTimes']

DataFileBlue_ = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataBlue-fitresults.npz')
paramsBlue_ = DataFileBlue_['FittedParams']
EclipseTimesBlue_ = DataFileBlue_['EclipseTimes']
DataFileBlue = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataBluefitresultsPaper.npz')
#DataFileBlue = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataBluefitresults.npz')
paramsBlue = DataFileBlue['FittedParams']
EclipseTimesBlue = DataFileBlue['EclipseTimes']
DataFileRed = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataRedfitresults.npz')
paramsRed = DataFileRed['FittedParams']
EclipseTimesRed = DataFileRed['EclipseTimes']

OldData = np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/CopperwheatResults.dat')
OldEclipseNum = np.array([line[0] for line in OldData])
OldEclipseTime = np.array([line[1] for line in OldData])
OldEclipseError = np.array([line[2] for line in OldData])

gap1=100
gap2=47
Eclipse=[0,1,2,gap1+2,gap1+3,gap1+4,gap1+4+gap2,gap1+5+gap2,gap1+6+gap2,gap1+7+gap2,gap1+8+gap2]
EclipseShort=[0,1,2,gap1+4+gap2,gap1+5+gap2,gap1+6+gap2,gap1+7+gap2,gap1+8+gap2] #doesn't include Dec10 eclipses

#EclipseTimes = [EclipseTimesAll,EclipseTimesBlue]
#EclipseNum = [Eclipse,EclipseShort]
#FixedEclipseTimes = [[],[]]
#FixedEclipseNum = [[],[]]

EclipseTimes = [EclipseTimesBlue]
EclipseNum = [EclipseShort]
FixedEclipseTimes = [[]]
FixedEclipseNum = [[]]

for ind in range(len(EclipseTimes)):
    for eclipse in range(len(EclipseTimes[ind])):
        if np.isnan(EclipseTimes[ind][eclipse])==False:
            FixedEclipseTimes[ind].append(EclipseTimes[ind][eclipse])
            FixedEclipseNum[ind].append(EclipseNum[ind][eclipse])

x0=125860
FixedEclipseNum = np.array([np.array(line)+x0 for line in FixedEclipseNum])

TotalEclipseNum = np.hstack((OldEclipseNum,FixedEclipseNum[0]))
TotalEclipseTime = np.hstack((OldEclipseTime,FixedEclipseTimes[0]))

OurEclipseNum = FixedEclipseNum[0]
OurEclipseTime = FixedEclipseTimes[0]

Std2006 = np.std(OldEclipseTime[:49]-linearEphem(OldEclipseNum[:49],ReportedOffset,ReportedPeriod))*24*60*60
Std2009 = np.std(OldEclipseTime[51:]-linearEphem(OldEclipseNum[51:],ReportedOffset,ReportedPeriod))*24*60*60

Std = np.std(OurEclipseTime-linearEphem(OurEclipseNum,ReportedOffset,ReportedPeriod))*24*60*60
ave = np.average(OurEclipseTime-linearEphem(OurEclipseNum,ReportedOffset,ReportedPeriod))*24*60*60

print 'std', Std2006
print Std2009,Std

sigma = np.hstack(([Std2006]*len(OldEclipseNum[:51]),np.hstack(([Std2009]*len(OldEclipseNum[51:]),[Std]*len(OurEclipseNum)))))


linearOpt = np.polyfit(TotalEclipseNum,TotalEclipseTime,1)
linearOffset = linearOpt[1]
linearPeriod = linearOpt[0]
print linearOpt

FixedOpt = fitFixedEphem(TotalEclipseTime,TotalEclipseNum,sigma)
QuadOffset = jd0 = FixedOpt[0]
f = FixedOpt[1]
#fDot = FixedOpt[2]
fDot=-1.2623365827*10**-9
QuadPdot = - fDot/(f**2)
QuadPeriod = 1/f
print "With fDot fixed:"
print f,fDot
#print 'offset %.15g'%QuadOffset
#print 'fixed period',QuadPeriod
print 'fixed pdot',QuadPdot

TimeDiff = []
for ind in range(len(TotalEclipseNum)):
    iNumber = TotalEclipseNum[ind]
    time = optimize.fsolve(fitModel,56000)
    TimeDiff.append((TotalEclipseTime[ind]-time[0])*24*3600)

ReducedChiSq3 = RedChiSq(TimeDiff,sigma,len(TotalEclipseTime)-3-1)
print 'red chi-sq',ReducedChiSq3

FixedOpt = fitEphem(TotalEclipseTime,TotalEclipseNum,sigma)
QuadOffset = jd0 = FixedOpt[0]
f = FixedOpt[1]
fDot = FixedOpt[2]
QuadPdot = - fDot/(f**2)
QuadPeriod = 1/f
print "With fDot fit using curvefit:"
print f,fDot
#print 'offset %.15g'%QuadOffset
#print 'fixed period',QuadPeriod
print 'fixed pdot',QuadPdot

TimeDiff = []
for ind in range(len(TotalEclipseNum)):
    iNumber = TotalEclipseNum[ind]
    time = optimize.fsolve(fitModel,56000)
    TimeDiff.append((TotalEclipseTime[ind]-time[0])*24*3600)

ReducedChiSq3 = RedChiSq(TimeDiff,sigma,len(TotalEclipseTime)-3-1)
print 'red chi-sq',ReducedChiSq3

    
coeff=np.polyfit(TotalEclipseTime,TotalEclipseNum,2)
a=coeff[2]
b=coeff[1]
c=coeff[0]
#N=f*(t-jd0)+fDot(t-jd0)**2
#N=(0.5*jd0**2*fDot-jd0*f)+(f-jd0*fDot)*t+(0.5*fDot)*t**2= a +b*t+ c*t**2
jd0=((b**2-4*a*c)**.5-b)/(2.*c)
fDot=2*c
f=(b**2-4*a*c)**.5
#print 'poly',jd0
QuadPeriod = 1/f
QuadPdot = - fDot/(f**2)
print "With fDot fit using polyfit:"
print f,fDot
print 'fixed pdot',QuadPdot

TimeDiff = []
for ind in range(len(TotalEclipseNum)):
    iNumber = TotalEclipseNum[ind]
    time = optimize.fsolve(fitModel,56000)
    TimeDiff.append((TotalEclipseTime[ind]-time[0])*24*3600)
ReducedChiSq3 = RedChiSq(TimeDiff,sigma,len(TotalEclipseTime)-3-1)
print 'red chi-sq',ReducedChiSq3

fig=plt.figure(figsize=(8.5,11))
plt.subplots_adjust(wspace=.3,hspace=.7)
font = {'size'   : 10}
matplotlib.rc('font', **font)


ReducedChiSq1 = RedChiSq((TotalEclipseTime-linearEphem(TotalEclipseNum,ReportedOffset,ReportedPeriod))*24*60*60,sigma,len(TotalEclipseTime)-2-1)

ax1=fig.add_subplot(331)
ax1.errorbar(OldEclipseNum[:49],(OldEclipseTime[:49]-linearEphem(OldEclipseNum[:49],ReportedOffset,ReportedPeriod))*24*60*60,yerr= [Std2006]*len(OldEclipseError[:49]),fmt='g^')
ax1.set_ylabel('O-C (s)')

ax2=fig.add_subplot(332, sharey=ax1)
ax2.errorbar(OldEclipseNum[51:],(OldEclipseTime[51:]-linearEphem(OldEclipseNum[51:],ReportedOffset,ReportedPeriod))*24*60*60,yerr= [Std2009]*len(OldEclipseError[51:]),fmt='g^')
ax2.set_title('Divergence from Reported Linear Ephemeris\n(Observed - Calculated)',fontsize=11)
ax2.xaxis.labelpad = 5
ax2.set_xlabel('Cycle Number')

ax3=fig.add_subplot(333, sharey=ax1)
ax3.errorbar(np.average(OurEclipseNum),ave,yerr=Std/(1.*np.sqrt(len(OurEclipseNum))),fmt='r*',label="average of our data = %f with std = %f\nReduced Chi-Squared (X^2)= %.3g"%(ave,Std,ReducedChiSq1))
ax3.errorbar(OurEclipseNum,(OurEclipseTime-linearEphem(OurEclipseNum,ReportedOffset,ReportedPeriod))*24*60*60,yerr=[Std]*len(OurEclipseNum),fmt='ro')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, borderaxespad=0.,prop={'size':'small'})
yticklabels = ax2.get_yticklabels()+ax3.get_yticklabels()
plt.setp(yticklabels, visible=False)



ReducedChiSq2 = RedChiSq((TotalEclipseTime-linearEphem(TotalEclipseNum,linearOffset,linearPeriod))*24*60*60,sigma,len(TotalEclipseTime)-2-1)

#Std2 = np.std(OurEclipseTime-Ephem(OurEclipseNum,linearPeriod,0))*24*60*60
ave2 = np.average(OurEclipseTime-linearEphem(OurEclipseNum,linearOffset,linearPeriod))*24*60*60

ax4=fig.add_subplot(334, sharex=ax1)
ax4.errorbar(OldEclipseNum[:49],(OldEclipseTime[:49]-linearEphem(OldEclipseNum[:49],linearOffset,linearPeriod))*24*60*60,yerr= [Std2006]*len(OldEclipseError[:49]),fmt='g^')
ax4.set_ylabel('O-C (s)')

ax5=fig.add_subplot(335, sharex=ax2,sharey=ax4)
ax5.errorbar(OldEclipseNum[51:],(OldEclipseTime[51:]-linearEphem(OldEclipseNum[51:],linearOffset,linearPeriod))*24*60*60,yerr= [Std2006]*len(OldEclipseError[51:]),fmt='g^')
plt.annotate('Divergence from Fitted Linear Ephemeris', xy=(0.5, 1.27), xycoords='axes fraction',ha='center',fontsize=11)
plt.annotate('Reported Period = 0.01966127289, Fitted Period = %.10f\nPeriod Difference (Reported-Fitted) = %.5g seconds\n t0 Difference = %.5g seconds '%(linearPeriod,(ReportedPeriod-linearPeriod)*24*60*60,(ReportedOffset-linearOffset)*24*60*60), xy=(0.5, 1.05), xycoords='axes fraction',ha='center',fontsize=9)
ax5.xaxis.labelpad = 5
ax5.set_xlabel('Cycle Number')

ax6=fig.add_subplot(336, sharex=ax3, sharey=ax4)
ax6.errorbar(np.average(OurEclipseNum),ave2,yerr=Std/(1.*np.sqrt(len(OurEclipseNum))),fmt='r*',label="average of our data = %f, X^2 = %.3g"%(ave2,ReducedChiSq2))
ax6.errorbar(OurEclipseNum,(OurEclipseTime-linearEphem(OurEclipseNum,linearOffset,linearPeriod))*24*60*60,yerr=[Std]*len(OurEclipseNum),fmt='ro')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, borderaxespad=0.,prop={'size':'small'})

yticklabels2 = ax5.get_yticklabels()+ax6.get_yticklabels()
plt.setp(yticklabels2, visible=False)



ReducedChiSq3 = RedChiSq(TimeDiff,sigma,len(TotalEclipseTime)-3-1)

#Std3 = np.std(TimeDiff[-len(OurEclipseNum):])
ave3 = np.average(TimeDiff[-len(OurEclipseNum):])

ax7=fig.add_subplot(337, sharex=ax1)
ax7.errorbar(OldEclipseNum[:49],TimeDiff[:49],yerr= [Std2006]*len(OldEclipseError[:49]),fmt='g^')
ax7.set_ylabel('O-C (s)')

ax8=fig.add_subplot(338, sharex=ax2,sharey=ax7)
ax8.errorbar(OldEclipseNum[51:],TimeDiff[51:-len(OurEclipseNum)],yerr= [Std2006]*len(OldEclipseError[51:]),fmt='g^')
plt.annotate('Divergence from Fitted Quadratic Ephemeris', xy=(0.5, 1.27), xycoords='axes fraction',ha='center',fontsize=11)
plt.annotate('Fitted Period = %.10f, Fitted Pdot = %.5g\nPeriod Difference (Reported-Fitted) = %.5g seconds\nt0 Difference = %.5g seconds'%(QuadPeriod,QuadPdot,(ReportedPeriod-QuadPeriod)*24*60*60,(ReportedOffset-QuadOffset)*24*60*60), xy=(0.5, 1.05), xycoords='axes fraction',ha='center',fontsize=9)
ax8.xaxis.labelpad = 5
ax8.set_xlabel('Cycle Number')

ax9=fig.add_subplot(339, sharex=ax3, sharey=ax7)
ax9.errorbar(np.average(OurEclipseNum),ave3,yerr=Std/(1.*np.sqrt(len(OurEclipseNum))),fmt='r*',label="average of our data = %f, X^2 = %.3g"%(ave3,ReducedChiSq3))
ax9.errorbar(OurEclipseNum,TimeDiff[-len(OurEclipseNum):],yerr=[Std]*len(OurEclipseNum),fmt='ro')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, borderaxespad=0.,prop={'size':'small'})

yticklabels3 = ax8.get_yticklabels()+ax9.get_yticklabels()
plt.setp(yticklabels3, visible=False)

plt.show()





