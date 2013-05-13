import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from scipy import pi

ReportedPeriod = 0.01966127289
#non barycentered JD: [2456271.0095138885081,2456271.0292013883591,2456271.0488194441423,2456273.0148495365866,2456273.0343981478363,
#2456273.0540740732104,2456273.9977083327249,2456274.0174305550754,2456274.0370833328925,2456274.0567824067548]

Dec8Eclipse1 = 56270.5136024671228 +3/2./3600/24 -17.584978980497891758816768/24./3600.
Dec8Eclipse2 = 56270.533291445470240 +3/2./3600/24 -17.584978980497891758816768/24./3600.
Dec8Eclipse3 = 56270.552910979538865 +3/2./3600/24 -17.584978980497891758816768/24./3600.
Dec10Eclipse1 = 56272.519043940577831 +10/2./3600/24 -17.584978980497891758816768/24./3600.
Dec10Eclipse2 = 56272.538617096048256 +10/2./3600/24 -17.584978980497891758816768/24./3600.
Dec10Eclipse3 = 56272.558294544789533 +10/2./3600/24 -17.584978980497891758816768/24./3600.
Dec11Eclipse0 = 56273.482310350409534 +10/2./3600/24 -17.584978980497891758816768/24./3600.
Dec11Eclipse1 = 56273.501987395095057 +10/2./3600/24 -17.584978980497891758816768/24./3600.
Dec11Eclipse2 = 56273.521699506782170 +10/2./3600/24 -17.584978980497891758816768/24./3600.
Dec11Eclipse3 = 56273.541376890119864 +10/2./3600/24 -17.584978980497891758816768/24./3600.
Dec11Eclipse4 = 56273.561088993898011 +10/2./3600/24 -17.584978980497891758816768/24./3600.
 
DataFileBlue = np.load('/Scratch/dataProcessing/SDSS_J0926/AllDataBluefitresults.npz')
paramsBlue = DataFileBlue['FittedParams']
EclipseTimesBlue = DataFileBlue['EclipseTimes']

gap1=100
gap2=47

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

def Ephem(x,Period,Pdot):
#    return 53795.9455191+Period*(x)+.5*(Pdot)*(x)**2

    Period = ReportedPeriod
    jd0 = 53795.9455191
    return jd0+Period*x/(1.-Pdot*x)


def fitEclipseNum(x,data):
    p0=(125860)
    popt, pcov = optimize.curve_fit(ephem, x, data, p0=p0)
    print p0
    return popt

def fitEphemeris(x,data):
    p0=(ReportedPeriod,-10**-7)
    popt, pcov = optimize.curve_fit(Ephem, x, data, p0=p0)
    print p0
    return popt

jd = np.array([Dec8Eclipse1,Dec8Eclipse2,Dec8Eclipse3,Dec10Eclipse1,Dec10Eclipse2,Dec10Eclipse3,Dec11Eclipse0,Dec11Eclipse1,Dec11Eclipse2,Dec11Eclipse3,Dec11Eclipse4])
Eclipse=[0,1,2,gap1+2,gap1+3,gap1+4,gap1+4+gap2,gap1+5+gap2,gap1+6+gap2,gap1+7+gap2,gap1+8+gap2]

MJDTDB = []

for eclipse in range(len(EclipseTimesBlue)):
    if EclipseTimesBlue[eclipse] != np.nan:
        MJDTDB.append(EclipseTimesBlue[eclipse])
    else:
        Eclipse.pop(eclipse)  

print Dec8Eclipse1,MJDTDB[0]
  
Eclipse = np.array(Eclipse)        
x0=125860
Eclipse=Eclipse+x0

#print Eclipse

#x0 = fitEclipseNum(Eclipse,jd)
#print x0
#x0=np.around(x0)

pOpt = fitEphemeris(Eclipse,MJDTDB)
print pOpt

fig=plt.figure()

ax=fig.add_subplot(311)
#ax.set_title(' Reported Period = 0.01966127289, Fitted Period = %.10f, Fitted Pdot = %.15f\nPeriod Difference (Reported-Fitted) = %.10f days or %.10f seconds\nAverage Period = %.10f, std of Period = %.10f days or %.10f seconds\nPeriod Difference (Reported - Average) = %.10f days or %.10f seconds'%(pOpt[0],pOpt[1],ReportedPeriod-pOpt[0],(ReportedPeriod-pOpt[0])*24*60*60,AveEph,std,std*24*60*60,ReportedPeriod-AveEph,(ReportedPeriod-AveEph)*24*3600))
#ax.set_title(' Reported Period = 0.01966127289\nAverage Period = %.10f, std of Period = %.10f days or %.10f seconds\nPeriod Difference (Reported - Average) = %.10f days or %.10f seconds'%(AveEph,std,std*24*60*60,ReportedPeriod-AveEph,(ReportedPeriod-AveEph)*24*3600))
#ax.plot(np.array(eclipseNum)+x0,times,'bo',label="data")
#ax.plot(np.array([3+x0,int(gap1)+2+x0]),[Gap108,Gap108],'k:')
#ax.plot(np.array([int(gap1)+5+x0,int(gap1)+4+int(gap2)+x0]),[Gap1110,Gap1110],'k:',label="Average Period Between Observations")
#ax.plot([1,int(gap1)+7+int(gap2)],[0.01966127289,0.01966127289],'r-')
#ax.plot([1+x0,int(gap1)+7+int(gap2)+x0],[AveEph,AveEph],'r*',label="Average Period")
#ax.set_ylabel('Period Between Eclipses')
#plt.legend(loc=4)

fullEclipseNum = np.arange(0,126000,1000.)
ax2=fig.add_subplot(312)
ax2.plot(Eclipse,MJDTDB,'ro')
ax2.plot(Eclipse,Ephem(Eclipse,ReportedPeriod,0),'b-',label="Ephemeris Formula")
ax2.plot(Eclipse,Ephem(Eclipse,pOpt[0],pOpt[1]),'r:',label="Ephemeris Fit")
ax2.set_ylabel('MJD(TDB)')
ax2.set_title('Data vs. Ephemeris Formula')

ax3=fig.add_subplot(313)
ax3.plot(Eclipse,(Ephem(Eclipse,ReportedPeriod,0)-MJDTDB)*24*60*60,'bo')
ax3.plot(Eclipse,(Ephem(Eclipse,ReportedPeriod,0)-jd)*24*60*60,'r.')
ax3.set_ylabel('O-C (s)')
ax3.set_title('Divergence from Reported Ephemeris: Equation - Data')

#ax4=fig.add_subplot(414)
#ax4.plot(Eclipse,(Ephem(Eclipse,pOpt[0],pOpt[1])-jd)*24*60*60,'ro')
plt.xlabel('Cycle')
#ax4.set_ylabel('O-C (s)')
#ax4.set_title('Divergence from Fitted Ephemeris')
plt.show()

def ephem(x,x0):
    return 53795.9455191+ReportedPeriod*(x+x0)
