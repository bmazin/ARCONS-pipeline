import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy import pi

ReportedPeriod = 0.01966127289
#non barycentered JD: [2456271.0095138885081,2456271.0292013883591,2456271.0488194441423,2456273.0148495365866,2456273.0343981478363,
#2456273.0540740732104,2456273.9977083327249,2456274.0174305550754,2456274.0370833328925,2456274.0567824067548]
#blue converted jd:
#np.array([56270.5136024671228,56270.533291445470240,56270.552910979538865,56272.519043940577831,56272.538617096048256, 56272.558294544789533,56273.482310350409534,56273.501987395095057,56273.521699506782170,56273.541376890119864, 56273.561088993898011]) +3/2./3600/24 -17.584978980497891758816768/24./3600.

OldData = np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/CopperwheatResults.dat')
OldEclipseNum = np.array([line[0] for line in OldData])
OldEclipseTime = np.array([line[1] for line in OldData])
OldEclipseError = np.array([line[2] for line in OldData])
 
DataFileBlue = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataBluefitresults.npz')
paramsBlue = DataFileBlue['FittedParams']
EclipseTimesBlue = DataFileBlue['EclipseTimes']
DataFileRed = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataRedfitresults.npz')
paramsRed = DataFileRed['FittedParams']
EclipseTimesRed = DataFileRed['EclipseTimes']
#DataFileMult = np.load('/Scratch/dataProcessing/SDSS_J0926/Dec8AddColorsfitresults.npz')
#paramsMult = DataFileMult['FittedParams']
#EclipseTimesMult = DataFileMult['EclipseTimes']
#EclipseNumMult = [0,1,2,0,1,2,0,1,2]

gap1=100
gap2=47


def Ephem(x,Period,Pdot):
#    return 53795.9455191+Period*(x)+.5*(Pdot)*(x)**2

#    Period = ReportedPeriod
    jd0 = 53795.9455191
    return jd0+Period*x/(1.-Pdot*x)

def fitEphemeris(x,data,sigma):
    p0=(ReportedPeriod,-10**-7)
    popt, pcov = optimize.curve_fit(Ephem, x, data, p0=p0,sigma=sigma)
    print p0
    return popt

jd = np.array([Dec8Eclipse1,Dec8Eclipse2,Dec8Eclipse3,Dec10Eclipse1,Dec10Eclipse2,Dec10Eclipse3,Dec11Eclipse0,Dec11Eclipse1,Dec11Eclipse2,Dec11Eclipse3,Dec11Eclipse4])
Eclipse=[0,1,2,gap1+2,gap1+3,gap1+4,gap1+4+gap2,gap1+5+gap2,gap1+6+gap2,gap1+7+gap2,gap1+8+gap2]

MJDTDB = []
OurEclipseNum = []
Dec8EclipseTimes = []
Dec8EclipseNum = []
Dec10EclipseTimes = []
Dec10EclipseNum = []
Dec11EclipseTimes = []
Dec11EclipseNum = []
EclipseTimes = [EclipseTimesBlue,EclipseTimesRed]
EclipseNum = [Eclipse,Eclipse]

for ind in range(len(EclipseTimes)):
    for eclipse in range(len(EclipseTimes[ind])):
        if np.isnan(EclipseTimes[ind][eclipse])==False:
#            print EclipseTimes[ind][eclipse]
            MJDTDB.append(EclipseTimes[ind][eclipse])
            OurEclipseNum.append(EclipseNum[ind][eclipse])
            if EclipseNum[ind][eclipse] < 3:
                Dec8EclipseTimes.append(EclipseTimes[ind][eclipse])
                Dec8EclipseNum.append(EclipseNum[ind][eclipse])
            elif EclipseNum[ind][eclipse] < 120:
                Dec10EclipseTimes.append(EclipseTimes[ind][eclipse])
                Dec10EclipseNum.append(EclipseNum[ind][eclipse])  
            else:
                Dec11EclipseTimes.append(EclipseTimes[ind][eclipse])
                Dec11EclipseNum.append(EclipseNum[ind][eclipse])

ii = np.where(np.array(OurEclipseNum) == 0)[0]
print [MJDTDB[index] for index in ii]
Dec81 = np.median([MJDTDB[index] for index in ii])
print Dec81
ii = np.where(np.array(OurEclipseNum) == 1)[0]
Dec82 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 2)[0]
Dec83 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 102)[0]
Dec101 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 103)[0]
Dec102 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 104)[0]
Dec103 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 151)[0]
Dec110 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 152)[0]
Dec111 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 153)[0]
Dec112 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 154)[0]
Dec113 = np.median([MJDTDB[index] for index in ii])
ii = np.where(np.array(OurEclipseNum) == 155)[0]
Dec114 = np.median([MJDTDB[index] for index in ii])

MedianEclipseTime = np.array([Dec81,Dec82,Dec83,Dec101,Dec102,Dec103,Dec110,Dec111,Dec112,Dec113,Dec114])

x0=125860
Eclipse = np.array(OurEclipseNum)+x0 
MJDTDB = np.array(MJDTDB)
Dec8EclipseTimes = np.array(Dec8EclipseTimes)
Dec8EclipseNum = np.array(Dec8EclipseNum)+x0
Dec10EclipseTimes = np.array(Dec10EclipseTimes)
Dec10EclipseNum = np.array(Dec10EclipseNum)+x0
Dec11EclipseTimes = np.array(Dec11EclipseTimes)
Dec11EclipseNum = np.array(Dec11EclipseNum)+x0     

#x0 = fitEclipseNum(Eclipse,jd)
#print x0
#x0=np.around(x0)
Std = np.std(MJDTDB-Ephem(Eclipse,ReportedPeriod,0))*24*60*60
ave = np.average((MJDTDB-Ephem(Eclipse,ReportedPeriod,0))*24*60*60)
Std8 = np.std(Dec8EclipseTimes-Ephem(Dec8EclipseNum,ReportedPeriod,0))*24*60*60
ave8 = np.average(Dec8EclipseTimes-Ephem(Dec8EclipseNum,ReportedPeriod,0))*24*60*60
Std10 = np.std(Dec10EclipseTimes-Ephem(Dec10EclipseNum,ReportedPeriod,0))*24*60*60
ave10 = np.average(Dec10EclipseTimes-Ephem(Dec10EclipseNum,ReportedPeriod,0))*24*60*60
Std11 = np.std(Dec11EclipseTimes-Ephem(Dec11EclipseNum,ReportedPeriod,0))*24*60*60
ave11 = np.average(Dec11EclipseTimes-Ephem(Dec11EclipseNum,ReportedPeriod,0))*24*60*60

print Std,ave
print Std8,ave8
print Std10,ave10
print Std11,ave11

TotalEclipseNum = np.hstack((np.hstack((OldEclipseNum,Dec8EclipseNum)),np.hstack((Dec10EclipseNum,Dec11EclipseNum))))

TotalEclipseTime = np.hstack((np.hstack((OldEclipseTime,Dec8EclipseTimes)),np.hstack((Dec10EclipseTimes,Dec11EclipseTimes))))

TotalEclipseError = np.hstack((OldEclipseError,np.array([Std8]*len(Dec8EclipseNum)+[Std10]*len(Dec10EclipseNum)+[Std11]*len(Dec11EclipseNum))))
print 'Eclipse Num = ',TotalEclipseNum
print 'Eclipse Times = ',TotalEclipseTime

pOpt = fitEphemeris(TotalEclipseNum,TotalEclipseTime,TotalEclipseError)
print pOpt

#StdFit = np.std(MJDTDB-Ephem(Eclipse,pOpt[0],pOpt[1]))*24*60*60
#aveFit = np.average((MJDTDB-Ephem(Eclipse,pOpt[0],pOpt[1]))*24*60*60)

fig=plt.figure()

#ax=fig.add_subplot(311)

#ax.set_title(' Reported Period = 0.01966127289\nAverage Period = %.10f, std of Period = %.10f days or %.10f seconds\nPeriod Difference (Reported - Average) = %.10f days or %.10f seconds'%(AveEph,std,std*24*60*60,ReportedPeriod-AveEph,(ReportedPeriod-AveEph)*24*3600))
#ax.plot(np.array(eclipseNum)+x0,times,'bo',label="data")
#ax.plot(np.array([3+x0,int(gap1)+2+x0]),[Gap108,Gap108],'k:')
#ax.plot(np.array([int(gap1)+5+x0,int(gap1)+4+int(gap2)+x0]),[Gap1110,Gap1110],'k:',label="Average Period Between Observations")
#ax.plot([1,int(gap1)+7+int(gap2)],[0.01966127289,0.01966127289],'r-')
#ax.plot([1+x0,int(gap1)+7+int(gap2)+x0],[AveEph,AveEph],'r*',label="Average Period")
#ax.set_ylabel('Period Between Eclipses')
#plt.legend(loc=4)

#fullEclipseNum = np.arange(0,126000,1000.)
#ax2=fig.add_subplot(211)
#ax2.plot(Eclipse,MJDTDB,'ro')
#ax2.plot(Eclipse,Ephem(Eclipse,ReportedPeriod,0),'b-',label="Ephemeris Formula")
#ax2.plot(Eclipse,Ephem(Eclipse,pOpt[0],pOpt[1]),'r:',label="Ephemeris Fit")
#ax2.set_ylabel('MJD(TDB)')
#ax2.set_title('Data vs. Ephemeris Formula')

ax3=fig.add_subplot(211)
ax3.set_title(' Reported Period = 0.01966127289, Fitted Period = %.10f, Fitted Pdot = %.15f\nPeriod Difference (Reported-Fitted) = %.10f days or %.10f seconds\nFit uses std as the error of our data\n\nDivergence from Reported Ephemeris: Data - Equation'%(pOpt[0],pOpt[1],ReportedPeriod-pOpt[0],(ReportedPeriod-pOpt[0])*24*60*60))
ax3.errorbar(OldEclipseNum,(OldEclipseTime-Ephem(OldEclipseNum,ReportedPeriod,0))*24*60*60, yerr=OldEclipseError, fmt='g^',label="Copperwheat Data")
ax3.plot(Eclipse,(MJDTDB-Ephem(Eclipse,ReportedPeriod,0))*24*60*60,'bo',label="Our Data (from 300-500nm)")
ax3.plot(Eclipse[:11],(jd-Ephem(Eclipse[:11],ReportedPeriod,0))*24*60*60,'r.',label = "Our Old Data")
ax3.errorbar(np.average(Eclipse),ave,yerr=Std,fmt='b>',label = "Average Difference of Our Data with total error = std = %f seconds"%(Std))
ax3.errorbar(np.average(Dec8EclipseNum),ave8,yerr=Std8,fmt='r>',label = "Average Difference of Our Dec 8 Data with error = std = %f seconds"%(Std8))
ax3.errorbar(np.average(Dec10EclipseNum),ave10,yerr=Std10,fmt='m>',label = "Average Difference of Our Dec 10 Data with error = std = %f seconds"%(Std10))
ax3.errorbar(np.average(Dec11EclipseNum),ave11,yerr=Std11,fmt='c>',label = "Average Difference of Our Dec 11 Data with error = std = %f seconds"%(Std11))
ax3.errorbar(np.average(Eclipse[:11]),np.average(MedianEclipseTime-Ephem(Eclipse[:11],ReportedPeriod,0))*24*3600,yerr=np.std(MedianEclipseTime-Ephem(Eclipse[:11],ReportedPeriod,0))*24*3600,fmt='ro',label = "Average Difference of Median Data with error = std = %f seconds"%(np.std(MedianEclipseTime-Ephem(Eclipse[:11],ReportedPeriod,0))*24*3600))
ax3.plot(Eclipse[:11],(MedianEclipseTime-Ephem(Eclipse[:11],ReportedPeriod,0))*24*60*60,'go')
ax3.plot(Eclipse[:11],(MedianEclipseTime-Ephem(Eclipse[:11],ReportedPeriod,0))*24*60*60,'go')
ax3.set_ylabel('O-C (s)')
#ax3.set_title('Divergence from Reported Ephemeris: Data - Equation, std = %f'%(np.std((jd-Ephem(Eclipse,ReportedPeriod,0))*24*60*60)))
plt.legend(loc='upper left',prop={'size':'small'})
ax4=fig.add_subplot(212)
ax4.plot(Eclipse,(MJDTDB-Ephem(Eclipse,pOpt[0],pOpt[1]))*24*60*60,'bo',label="Our Data (from 300-500nm)")
ax4.plot(OldEclipseNum,(OldEclipseTime-Ephem(OldEclipseNum,pOpt[0],pOpt[1]))*24*60*60,'g^',label="Copperwheat Data")
#ax4.errorbar(np.average(Eclipse),aveFit,yerr=StdFit,fmt='b>',label ="Average Difference of Our Data with error = std = %f seconds"%(StdFit))
plt.xlabel('Cycle')
ax4.set_ylabel('O-C (s)')
ax4.set_title('Divergence from Fitted Ephemeris')
plt.legend(loc='upper left',prop={'size':'small'})
plt.show()

##############################


def fitEclipseNum(x,data):
    p0=(125860)
    popt, pcov = optimize.curve_fit(ephem, x, data, p0=p0)
    print p0
    return popt

#MJD =JD-2400000.5
#MJDtdb = MJDtcb+(-1.550519768*10**-8*(MJDtcb+2400000.5-2443144.5003725)*86400-6.55*10**-5)/(24.*3600.)
#(56270.5136024671228-53795.9455191)/0.01966127289=125860.014110567731
#fit and x0=125860.01136206
#times = [Dec8Eclipse2-Dec8Eclipse1,Dec8Eclipse3-Dec8Eclipse2,Dec10Eclipse2-Dec10Eclipse1,Dec10Eclipse3-Dec10Eclipse2,Dec11Eclipse1-Dec11Eclipse0,Dec11Eclipse2-Dec11Eclipse1,Dec11Eclipse3-Dec11Eclipse2,Dec11Eclipse4-Dec11Eclipse3]
#print times
#eph = np.average(times)
#std = np.std(times)
#print 'look', np.average(Dec8Eclipse2-Dec8Eclipse1,Dec8Eclipse3-Dec8Eclipse2)
#print eph,std,std*24*60*60
#gap1 = (Dec10Eclipse1-Dec8Eclipse3)/eph
#gap2 = (Dec11Eclipse0-Dec10Eclipse3)/eph
#print gap1,gap2
#gap1=np.around(gap1)
#gap2=np.around(gap2)
#Gap108 = (Dec10Eclipse1-Dec8Eclipse3)/gap1
#Gap1110= (Dec11Eclipse0-Dec10Eclipse3)/gap2
#print Gap108, Gap1110
#AveEph = np.average(times+[Gap108]*gap1+[Gap1110]*gap2) 
#print AveEph
#eclipseNum = [1,2,gap1+3,gap1+4,gap1+5+gap2,gap1+6+gap2,gap1+7+gap2,gap1+8+gap2]
#print eclipseNum,times

def ephem(x,x0):
    return 53795.9455191+ReportedPeriod*(x+x0)
