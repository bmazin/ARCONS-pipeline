import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy import pi
import matplotlib

###This program goes through a range of P-dot values and finds the quality of the fit (reduced chi-squared) for each one. It then plots these points with a line indicating what P-dots are within one sigma of the reduced chi-squared.

### At the end the commented out region is used to find the max and min P-dot values within one sigma of the reduced chi-squared. This currently deletes data and gives an error when it tries to plot. A crappy work around is to run it once not commented out to get the limits, and a second time with it commented out to get the plot.

def FixedEphem(t,jd0,f,fDot):
    return f*(t-jd0) + 0.5 * fDot * (t-jd0)**2

def fitModel(t):
    return iNumber - FixedEphem(t,jd0,f,fDot)

def fitFixedEphem(t,data,sigma):
    p0=(ReportedOffset,1/ReportedPeriod,-10**-13/ReportedPeriod**2)
    popt, pcov = optimize.curve_fit(FixedEphem, t, data, p0=p0,sigma=sigma)
#    print 'quad stuff', pcov
#    print p0
    return popt

def RedChiSq(diff,sigma,DoF):
    return np.sum([diff[ind]**2/(1.*sigma[ind]**2*DoF) for ind in range(len(diff))]) 
    #DoF = N(num of points)-v(num variables)-1

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
EclipseShort=[0,1,2,gap1+4+gap2,gap1+5+gap2,gap1+6+gap2,gap1+7+gap2,gap1+8+gap2]

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


FixedOpt = fitFixedEphem(TotalEclipseTime,TotalEclipseNum,sigma)
QuadOffset = jd0 = FixedOpt[0]
f = FixedOpt[1]
fDot = FixedOpt[2]
QuadPdot = - fDot/(f**2)
QuadPeriod = 1/f
print 'offset %.15g'%QuadOffset
print 'fixed period',QuadPeriod
print 'fixed pdot',QuadPdot

TimeDiff2=[]
for ind in range(len(TotalEclipseNum)):
    iNumber = TotalEclipseNum[ind]
    time = optimize.fsolve(fitModel,56000)
    TimeDiff2.append((TotalEclipseTime[ind]-time[0])*24*3600)

ReducedChiSq = RedChiSq(TimeDiff2,sigma,len(TotalEclipseTime)-3-1)

print 'chisq',ReducedChiSq

Pdots = []
ChiSqs = []
Ps = []
Offsets = []

def fitFixedFdotEphem(t,data,sigma):
    p0=(ReportedOffset,1/ReportedPeriod)
    popt, pcov = optimize.curve_fit(FixedFdotEphem, t, data, p0=p0,sigma=sigma)
    return popt

def fitFixedFdotModel(t):
    return iNumber - FixedFdotEphem(t,jd0,f)

for SetPDot in np.array(np.linspace(-15,25,500))*10**(-13):
    fDot = SetfDot = -SetPDot/(QuadPeriod)**2
    def FixedFdotEphem(t,jd0,f):
        fDot = SetfDot
        return f*(t-jd0) + 0.5 * fDot * (t-jd0)**2
    pOpt = fitFixedFdotEphem(TotalEclipseTime,TotalEclipseNum,sigma)
    
    QuadOffset = jd0 = pOpt[0]
    f = pOpt[1]
    QuadPdot = - fDot/(f**2)
    QuadPeriod = 1/f

    TimeDiff = []
    for ind in range(len(TotalEclipseNum)):
        iNumber = TotalEclipseNum[ind]
        time = optimize.fsolve(fitFixedFdotModel,56000)
        TimeDiff.append((TotalEclipseTime[ind]-time[0])*24*3600)

    ReducedChiSq3 = RedChiSq(TimeDiff,sigma,len(TotalEclipseTime)-3-1)
    Pdots.append(QuadPdot)
    ChiSqs.append(ReducedChiSq3)
    Ps.append(QuadPeriod)
    Offsets.append(QuadOffset)

indexmin=np.array(ChiSqs).argmin()
minChiSq = ChiSqs[indexmin]
minPdot = Pdots[indexmin]
print 'minimum reduced chi-squared',minChiSq
print 'optimum p-dot found by minimizing red chi-sq', minPdot
print Ps[indexmin],Offsets[indexmin]

print 'Number of Points',len(TotalEclipseTime)

ChiError = np.sqrt(2/(1.*len(TotalEclipseTime)))
print 'chi error %.5g'%ChiError

### This commented out region is used to find the max and min P-dot values within one sigma of the reduced chi-squared. This currently deletes data and gives an error when it tries to plot. A crappy work around is to run it once not commented out to get the limits, and a second time with it commented out to get the plot.


#ChiSqs2 = ChiSqs

#pdot = -1
#while pdot < 0:
#    idx = (np.abs(np.array(ChiSqs)-(minChiSq+ChiError))).argmin()
#    pdot = Pdots[idx]
#    print pdot
#    ChiSqs.pop(idx)
#    print len(ChiSqs)

#print len(Pdots),len(ChiSqs)
#print len(ChiSqs2)

fig = plt.figure()
plt.plot(Pdots,ChiSqs,'r-')
plt.plot([Pdots[0],Pdots[len(Pdots)-1]],[minChiSq+ChiError]*2,'b:')
plt.xlabel('Pdot')
plt.ylabel('Reduced Chi-Squared')
plt.minorticks_on
plt.show()





