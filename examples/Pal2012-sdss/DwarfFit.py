import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from util import utils
import scipy.optimize as optimize
from scipy import pi,sqrt,exp


#### This program fits, plots, and saves the fits for multiple eclipses simultaneously. For the sdss data, each eclipse's paramaters are shown for All wvl, Blue, and Red. Change what's included in files (line 79) to control what eclipses are fit and what wvl band. fileNames are used to title the plots. plotdimy and plotdimx (line 88) control how many subplots (one for each eclipse) and their orientation. The save is currently commented out at the end of the program, so you don't save over the last npz file saved. The npz file contains EclipseTimes and FittedParams. If the eclipse is not fit, nan is saved.


#WDp0=(1100,b0,np.average(data),(np.min(data)-np.max(data))/2.,1)

#[npzDataFile,Int,BMJDfile,fitwidth,RatioBelow,RatioAbove,EclipseNum, cutoffmin,cutoffmax,b0,peakCutOffmin,peakCutOffmax,fit?]

Dec8txtBjd = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3MJD_TDB.txt'

Dec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAllInt3newwvlcal.npz'
Dec8Eclipse1 = [Dec8npzData,3,Dec8txtBjd,23,1/5.,4/5.,1,0,300,0,0,300,True]
Dec8Eclipse2 = [Dec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,200,True] 
Dec8Eclipse3 = [Dec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,0,300,True]

BDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3.npz'
BDec8Eclipse1 = [BDec8npzData,3,Dec8txtBjd,23,1/5.,4/5.,1,0,300,0,0,300,True]
BDec8Eclipse2 = [BDec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,200,True] 
BDec8Eclipse3 = [BDec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,0,300,True]

#Int=10. #so 23 => 7, 300 => 90
#BDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10.npz'
Dec10txtBjd = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10MJD_TDB.txt'

Dec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfAllInt10newwvlcal.npz'
Dec10Eclipse1 = [Dec10npzData,10,Dec10txtBjd,7,1/3.,2/3.,1,490,580,.5,0,90,False]
Dec10Eclipse2 = [Dec10npzData,10,Dec10txtBjd,7,3/7.,2/3.,2,650,740,-.5,0,80,False]
Dec10Eclipse3 = [Dec10npzData,10,Dec10txtBjd,7,2/5.,2/3.,3,820,910,.3,0,90,False]


BDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlue-Int10newwvlcal.npz'
BDec10Eclipse1 = [BDec10npzData,10,Dec10txtBjd,7,1/3.,2/3.,1,490,580,.5,0,90,False]
BDec10Eclipse2 = [BDec10npzData,10,Dec10txtBjd,7,2/5.,2/3.,2,650,740,-.5,0,80,True]
BDec10Eclipse3 = [BDec10npzData,10,Dec10txtBjd,7,2/5.,2/3.,3,820,910,.3,0,90,True]


Dec11txtBjd = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10MJD_TDB.txt'

Dec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfAllInt10newwvlcal.npz'
Dec11Eclipse0 = [Dec11npzData,10,Dec11txtBjd,6,0.,1/5.,0,0,90,5.5,0,40,True] 
Dec11Eclipse1 = [Dec11npzData,10,Dec11txtBjd,7,3/7.,3/5.,1,130,220,.5,0,90,True]
Dec11Eclipse2 = [Dec11npzData,10,Dec11txtBjd,7,3/7.,4/7.,2,300,390,.5,0,90,True]
Dec11Eclipse3 = [Dec11npzData,10,Dec11txtBjd,7,2/5.,3/5.,3,470,560,.7,20,90,True]
Dec11Eclipse4 = [Dec11npzData,10,Dec11txtBjd,6,2/5.,3/5.,4,640,730,.7,40,90,False]

BDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10.npz'
BDec11Eclipse0 = [BDec11npzData,10,Dec11txtBjd,6,0.,1/5.,0,0,90,5.5,0,40,True] 
BDec11Eclipse1 = [BDec11npzData,10,Dec11txtBjd,7,3/7.,3/5.,1,130,220,.5,0,90,True]
BDec11Eclipse2 = [BDec11npzData,10,Dec11txtBjd,7,1/3.,2/3.,2,300,390,.5,0,90,True]
BDec11Eclipse3 = [BDec11npzData,10,Dec11txtBjd,7,2/5.,3/5.,3,470,560,.7,20,90,True]
BDec11Eclipse4 = [BDec11npzData,10,Dec11txtBjd,6,2/5.,3/5.,4,640,730,.7,40,90,False]

RDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfRedInt3.npz'
RDec8Eclipse1 = [RDec8npzData,3,Dec8txtBjd,23,1/5.,4/5.,1,3,303,0,0,300,True]
RDec8Eclipse2 = [RDec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,200,True] 
RDec8Eclipse3 = [RDec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,0,300,True]

RDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfRedInt10.npz'
RDec10Eclipse1 = [RDec10npzData,10,Dec10txtBjd,7,1/4.,3/4.,1,490,580,.5,40,70,True]
RDec10Eclipse2 = [RDec10npzData,10,Dec10txtBjd,7,3/8.,2/3.,2,650,740,-.5,25,80,False]
RDec10Eclipse3 = [RDec10npzData,10,Dec10txtBjd,7,2/5.,2/3.,3,820,910,.3,0,90,True]

RDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfRedInt10.npz'
RDec11Eclipse0 = [RDec11npzData,10,Dec11txtBjd,6,0.,1/3.,0,0,90,5.5,0,40,True] 
RDec11Eclipse1 = [RDec11npzData,10,Dec11txtBjd,7,2/5.,3/5.,1,130,220,.5,0,90,True]
RDec11Eclipse2 = [RDec11npzData,10,Dec11txtBjd,7,2/5.,3/5.,2,300,390,.5,0,90,True]
RDec11Eclipse3 = [RDec11npzData,10,Dec11txtBjd,7,2/5.,3/5.,3,470,560,.7,0,90,False]
RDec11Eclipse4 = [RDec11npzData,10,Dec11txtBjd,7,2/5.,3/5.,4,640,730,.7,0,90,False]



#files = [Dec8Eclipse1]
#files = [Dec8Eclipse1,Dec8Eclipse2,Dec8Eclipse3,Dec10Eclipse1,Dec10Eclipse2, Dec10Eclipse3,Dec11Eclipse0,Dec11Eclipse1,Dec11Eclipse2, Dec11Eclipse3,Dec11Eclipse4]
fileNames = ['Dec8Eclipse1','Dec8Eclipse2','Dec8Eclipse3','Dec11Eclipse0','Dec11Eclipse1','Dec11Eclipse2', 'Dec11Eclipse3','Dec11Eclipse4']#'Dec10Eclipse1','Dec10Eclipse2', 'Dec10Eclipse3',
files = [BDec8Eclipse1,BDec8Eclipse2,BDec8Eclipse3,BDec11Eclipse0,BDec11Eclipse1,BDec11Eclipse2,BDec11Eclipse3, BDec11Eclipse4]
#BDec10Eclipse1,BDec10Eclipse2,BDec10Eclipse3,
#files = [RDec8Eclipse1,RDec8Eclipse2,RDec8Eclipse3,RDec10Eclipse1,RDec10Eclipse2,RDec10Eclipse3,RDec11Eclipse0,RDec11Eclipse1,RDec11Eclipse2,RDec11Eclipse3,RDec11Eclipse4]
#files = [GDec8Eclipse1,GDec8Eclipse2,GDec8Eclipse3,R_Dec8Eclipse1,R_Dec8Eclipse2,R_Dec8Eclipse3,I_Dec8Eclipse1,I_Dec8Eclipse2,I_Dec8Eclipse3]

#files = [IDec8Eclipse1,IDec8Eclipse2,IDec8Eclipse3,IDec10Eclipse1,IDec10Eclipse2,IDec10Eclipse3,IDec11Eclipse0,IDec11Eclipse1,IDec11Eclipse2,IDec11Eclipse3,IDec11Eclipse4]
plotdimy = 3
plotdimx = 3
#plotdimy = 1
#plotdimx = 1

def zero_if_negative(x):
    if isinstance(x, (int, long, float, complex))==True:
        if x < 0:
            x=0
    else: 
        for item in range(len(x)):
             if x[item] < 0 or item <len(x)*RatioBelow or item > len(x)*RatioAbove:
                 x[item]=0
    return x

def removeOutliers(jd,data):
    ave = np.average(data)
    newjd=[]
    newdata=[]
    for index in range(len(data)):
        if data[index]> ave-400 and data[index] < ave+400 and data[index] > 5:
            newjd.append(jd[index])
            newdata.append(data[index])
    jd=np.array(newjd)
    data=np.array(newdata)
    return jd,data

#theta =< pi/2 which implies mu = cos(theta) => 0
#epsilon is a limb-darkening coefficient
#trapezoidally weighted averaging?
#try singular-value decomposition (SVD)= all components contribute linearly to the light curve
def WhiteDwarf(x,a,b,c,amp,epsilon):
    theta=a*x+b
    mu=np.cos(theta)
    mu=zero_if_negative(mu)
    return amp*mu*(1-epsilon+epsilon*mu)+c # Intensity is proportional to mu*(1-epsilon+epsilon*mu)

def fitWhiteDwarf(x,data): 
    sigma= [10]*fitstart+[1]*(fitend-fitstart)+[10]*(len(data)-fitend)
    popt, pcov = optimize.curve_fit(WhiteDwarf, x, data, p0=WDp0,sigma=sigma)
#    print WDp0
#    print pcov
    return popt

fig = plt.figure()#figsize=(8.5,11))

font = {'size'   : 8}
matplotlib.rc('font', **font)
plt.subplots_adjust(hspace=.5,wspace=.6)


EclipseTimes = []
FittedParams = []

for eclipse in range(len(files)):
    npzDataFile = files[eclipse][0]
    Int = files[eclipse][1]
    BMJDfile = files[eclipse][2]
    fitwidth = files[eclipse][3]
    RatioBelow = files[eclipse][4]
    RatioAbove = files[eclipse][5]
    EclipseNum = files[eclipse][6] 
    cutoffmin = files[eclipse][7]
    cutoffmax = files[eclipse][8]
    b0 = files[eclipse][9]
    peakCutOffmin = files[eclipse][10]
    peakCutOffmax = files[eclipse][11]
    fit = files[eclipse][12]

    BMJD = np.loadtxt(BMJDfile)
    DataFile = np.load(npzDataFile)
    params = DataFile['params']
#    jd = DataFile['jd']
    amps = params[:,1]
    widths = params[:,4]
    xpos = params[:,2]
    ypos = params[:,3]
    curve = 2*pi*amps*widths**2/Int
    curve /= np.average(curve)
    print len(curve)
    jd = BMJD
    jd = jd[cutoffmin:cutoffmax]
    data=curve[cutoffmin:cutoffmax]

    jd,data=removeOutliers(jd,data)
    x=jd-jd[0]

    WDp0=(1100,b0,np.average(data),(np.min(data)-np.max(data))/2.,1)

    DataIndex_min = data[peakCutOffmin:peakCutOffmax].argmin()+peakCutOffmin
    fitstart = DataIndex_min-fitwidth
    fitend = fitstart+2*fitwidth
    if fitend > len(jd)-1:
        fitend = len(jd)-1

    ax = fig.add_subplot(plotdimy,plotdimx,eclipse+1)
#    ax.xaxis.labelpad = 20
#    ax.yaxis.labelpad = 20

    print fitstart,fitend
## If fit fails, uncomment the next two comments so the suggested fit and the data are shown.

#    ax.plot(jd,WhiteDwarf(x, WDp0[0],WDp0[1],WDp0[2],WDp0[3],WDp0[4]),'--')
    ax.plot(jd,data,'k')
#    plt.show()


    if fit == True:
        ax.plot([jd[fitstart],jd[fitend]],[np.average(data)]*2,'go')
        TFparams = fitWhiteDwarf(x,data)
        print TFparams

        FittedParams.append(TFparams)
        Xrange = np.arange(x[0],x[len(x)-1]+.01/(24*3600),.01/(24*3600))

        WD = []
        for index in range(len(Xrange)):
            if index <len(Xrange)*RatioBelow or index > len(Xrange)*RatioAbove:
                WD.append(TFparams[2])
            else:
                WD.append(WhiteDwarf(Xrange[index], TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4]))

        WD= np.array(WD)
        WDindex_min=WD.argmin()
        WDx_min=Xrange[WDindex_min]
        WDeclipseTime = WDx_min+jd[0] 
        print 'Eclipse Time is %.15f'%(WDeclipseTime)

        EclipseTimes.append(WDeclipseTime)

        ax.set_title(np.str(fileNames[eclipse])[:3]+' '+ fileNames[eclipse][3:-8]+ ' eclipse %d\n'%(EclipseNum)+ 'Fitted Eclipse Time = %.13f'%(WDeclipseTime))#,size='small')

        ax.plot(Xrange+jd[0],WD,'r')#, label = 'White Dwarf Light Curve Fitted',)
    else:
        FittedParams.append([np.nan])
        EclipseTimes.append(np.nan)
        ax.set_title(fileNames[eclipse][:3]+' '+ fileNames[eclipse][3:-8]+ ' eclipse %d\n'%(EclipseNum))#,size='small')
    ax.set_xlabel('MJD(TDB)')#,size='small')
    ax.set_ylabel('Normalized Flux')

#np.savez('/Scratch/dataProcessing/SDSS_J0926/AllDataBluefitresults.npz',EclipseTimes = EclipseTimes, FittedParams = FittedParams)

plt.show()


