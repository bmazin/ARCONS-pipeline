import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from util import utils
import scipy.optimize as optimize
from scipy import pi,sqrt,exp

#WDp0=(1100,b0,np.average(data),(np.min(data)-np.max(data))/2.,1)

#[npzDataFile,Int,BMJDfile,fitwidth,RatioBelow,RatioAbove,EclipseNum, cutoffmin,cutoffmax,b0,peakCutOffmin,peakCutOffmax,fit?]
#BDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3.npz'
#BDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3.npz'
BDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlue-Int3newwvlcal.npz'
Dec8txtBjd = BDec8txtBjd = RDec8txtBjd = IDec8txtBjd = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3MJD_TDB.txt'
BDec8Eclipse1 = [BDec8npzData,3,BDec8txtBjd,23,1/5.,4/5.,1,0,300,0,0,300,True]
BDec8Eclipse2 = [BDec8npzData,3,BDec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,200,True] 
BDec8Eclipse3 = [BDec8npzData,3,BDec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,0,300,True]

#Int=10. #so 23 => 7, 300 => 90
#BDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10.npz'
BDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlue-Int10newwvlcal.npz'
BDec10txtBjd = RDec10txtBjd = IDec10txtBjd = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfBlueInt10MJD_TDB.txt'
BDec10Eclipse1 = [BDec10npzData,10,BDec10txtBjd,7,1/3.,2/3.,1,490,580,.5,0,90,False]
BDec10Eclipse2 = [BDec10npzData,10,BDec10txtBjd,7,2/5.,2/3.,2,650,740,-.5,0,80,True]
BDec10Eclipse3 = [BDec10npzData,10,BDec10txtBjd,7,2/5.,2/3.,3,820,910,.3,0,90,True]

#BDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10.npz'
BDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlue-Int10newwvlcal.npz'
BDec11txtBjd = RDec11txtBjd = IDec11txtBjd = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfBlueInt10MJD_TDB.txt'
BDec11Eclipse0 = [BDec11npzData,10,BDec11txtBjd,6,0.,1/5.,0,0,90,5.5,0,40,True] 
BDec11Eclipse1 = [BDec11npzData,10,BDec11txtBjd,7,3/7.,3/5.,1,130,220,.5,0,90,True]
BDec11Eclipse2 = [BDec11npzData,10,BDec11txtBjd,7,1/3.,2/3.,2,300,390,.5,0,90,False]
BDec11Eclipse3 = [BDec11npzData,10,BDec11txtBjd,7,2/5.,3/5.,3,470,560,.7,20,90,True]
BDec11Eclipse4 = [BDec11npzData,10,BDec11txtBjd,6,2/5.,3/5.,4,640,730,.7,40,90,False]

RDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfRedInt3.npz'
#RDec8txtBjd = '/Scratch/dataProcessing/SDSS_J0926/OlderData/Dec8fitpsfBlueMJD-TDB.txt'
RDec8Eclipse1 = [RDec8npzData,3,RDec8txtBjd,23,1/5.,4/5.,1,3,303,0,0,300,True]
RDec8Eclipse2 = [RDec8npzData,3,RDec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,200,True] 
RDec8Eclipse3 = [RDec8npzData,3,RDec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,0,300,True]

RDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfRedInt10.npz'
#RDec10txtBjd = '/Scratch/dataProcessing/SDSS_J0926/OlderData/Dec10fitpsfBlueInt10MJD-TDB.txt'
RDec10Eclipse1 = [RDec10npzData,10,RDec10txtBjd,7,1/4.,3/4.,1,490,580,.5,40,70,True]
RDec10Eclipse2 = [RDec10npzData,10,RDec10txtBjd,7,3/8.,2/3.,2,650,740,-.5,25,80,False]
RDec10Eclipse3 = [RDec10npzData,10,RDec10txtBjd,7,2/5.,2/3.,3,820,910,.3,0,90,True]

RDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfRedInt10.npz'
#RDec11txtBjd = '/Scratch/dataProcessing/SDSS_J0926/OlderData/Dec11fitpsfBlueInt10MJD-TDB.txt'
RDec11Eclipse0 = [RDec11npzData,10,RDec11txtBjd,6,0.,1/3.,0,0,90,5.5,0,40,True] 
RDec11Eclipse1 = [RDec11npzData,10,RDec11txtBjd,7,2/5.,3/5.,1,130,220,.5,0,90,True]
RDec11Eclipse2 = [RDec11npzData,10,RDec11txtBjd,7,2/5.,3/5.,2,300,390,.5,0,90,True]
RDec11Eclipse3 = [RDec11npzData,10,RDec11txtBjd,7,2/5.,3/5.,3,470,560,.7,0,90,False]
RDec11Eclipse4 = [RDec11npzData,10,RDec11txtBjd,7,2/5.,3/5.,4,640,730,.7,0,90,False]

GDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfGreenInt3.npz'
GDec8Eclipse1 = [GDec8npzData,3,Dec8txtBjd,23,1/5.,4/5.,1,3,303,0,0,300,True]
GDec8Eclipse2 = [GDec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,200,True] 
GDec8Eclipse3 = [GDec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,100,300,True]

R_Dec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfRed-Int3.npz'
R_Dec8Eclipse1 = [R_Dec8npzData,3,Dec8txtBjd,23,1/5.,4/5.,1,3,303,0,0,300,True]
R_Dec8Eclipse2 = [R_Dec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,180,True] 
R_Dec8Eclipse3 = [R_Dec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,100,300,True]

I_Dec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfInfra-Int3.npz'
I_Dec8Eclipse1 = [I_Dec8npzData,3,Dec8txtBjd,23,1/5.,4/5.,1,3,303,0,0,300,False]
I_Dec8Eclipse2 = [I_Dec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,200,False] 
I_Dec8Eclipse3 = [I_Dec8npzData,3,Dec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,100,300,False]

#IDec8npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfInfraInt3.npz'
#IDec8Eclipse1 = [IDec8npzData,3,IDec8txtBjd,23,1/5.,4/5.,1,3,303,0,0,300,False]
#IDec8Eclipse2 = [IDec8npzData,3,IDec8txtBjd,23,2/5.,4/5.,2,550,850,-.5,0,200,False] 
#IDec8Eclipse3 = [IDec8npzData,3,IDec8txtBjd,23,2/5.,4/5.,3,1100,1400,-1,100,300,False]

#IDec10npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec10fitpsfInfraInt10.npz'
#IDec10txtBjd = '/Scratch/dataProcessing/SDSS_J0926/OlderData/Dec10fitpsfBlueInt10MJD-TDB.txt'
#IDec10Eclipse1 = [IDec10npzData,10,IDec10txtBjd,7,1/4.,3/4.,1,490,580,.5,40,70,False]
#IDec10Eclipse2 = [IDec10npzData,10,IDec10txtBjd,9,1/3.,2/3.,2,650,740,-.5,25,80,False]
#IDec10Eclipse3 = [IDec10npzData,10,IDec10txtBjd,9,1/3.,2/3.,3,820,910,.3,0,90,False]

#IDec11npzData = '/Scratch/dataProcessing/SDSS_J0926/Dec11fitpsfInfraInt10.npz'
#IDec11txtBjd = '/Scratch/dataProcessing/SDSS_J0926/OlderData/Dec11fitpsfBlueInt10MJD-TDB.txt'
#IDec11Eclipse0 = [IDec11npzData,10,IDec11txtBjd,6,0.,1/3.,0,0,90,5.5,0,40,False] 
#IDec11Eclipse1 = [IDec11npzData,10,IDec11txtBjd,9,2/5.,3/5.,1,130,220,.5,0,90,False]
#IDec11Eclipse2 = [IDec11npzData,10,IDec11txtBjd,7,3/7.,3/5.,2,300,390,.5,0,90,False]
#IDec11Eclipse3 = [IDec11npzData,10,IDec11txtBjd,7,2/5.,3/5.,3,470,560,.7,0,90,False]
#IDec11Eclipse4 = [IDec11npzData,10,IDec11txtBjd,7,2/5.,3/5.,4,640,730,.7,0,90,False]

files = [BDec8Eclipse1,BDec8Eclipse2,BDec8Eclipse3,BDec10Eclipse1,BDec10Eclipse2,BDec10Eclipse3,BDec11Eclipse0,BDec11Eclipse1,BDec11Eclipse2,BDec11Eclipse3,BDec11Eclipse4]
#files = [RDec8Eclipse1,RDec8Eclipse2,RDec8Eclipse3,RDec10Eclipse1,RDec10Eclipse2,RDec10Eclipse3,RDec11Eclipse0,RDec11Eclipse1,RDec11Eclipse2,RDec11Eclipse3,RDec11Eclipse4]
#files = [GDec8Eclipse1,GDec8Eclipse2,GDec8Eclipse3,R_Dec8Eclipse1,R_Dec8Eclipse2,R_Dec8Eclipse3,I_Dec8Eclipse1,I_Dec8Eclipse2,I_Dec8Eclipse3]

#files = [IDec8Eclipse1,IDec8Eclipse2,IDec8Eclipse3,IDec10Eclipse1,IDec10Eclipse2,IDec10Eclipse3,IDec11Eclipse0,IDec11Eclipse1,IDec11Eclipse2,IDec11Eclipse3,IDec11Eclipse4]
plotdimy = 4
plotdimx = 3

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

fig = plt.figure()

font = {'size'   : 8}
matplotlib.rc('font', **font)
plt.subplots_adjust(hspace=.4,wspace=.6)

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
#jd = DataFile['jd']
#params[:,0] are the height offset of the gaussian (so the number of background pulses)
    amps = params[:,1]
    widths = params[:,4]
    xpos = params[:,2]
    ypos = params[:,3]
    curve = 2*pi*amps*widths**2/Int
    print len(curve)
#    jd = [line[1] for line in BMJD]
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
    print fitstart,fitend
    ax.plot(jd,WhiteDwarf(x, WDp0[0],WDp0[1],WDp0[2],WDp0[3],WDp0[4]),'--')
    ax.plot([jd[fitstart],jd[fitend]],[np.average(data)]*2,'go')
#    ax.plot([jd[peakCutOffmin],jd[peakCutOffmax-1]],[np.average(data)]*2,'ro')
#    ax.plot([jd[DataIndex_min]],[np.average(data)],'bo')
    ax.plot(jd,data,'k')
#    plt.show()

    if fit == True:
        print 'x',len(x)
        print 'data',len(data)
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
#        WDx = WhiteDwarf(x,TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4])

        WDindex_min=WD.argmin()
        WDx_min=Xrange[WDindex_min]
        WDeclipseTime = WDx_min+jd[0] 
        print 'Eclipse Time is %.15f'%(WDeclipseTime)

        EclipseTimes.append(WDeclipseTime)

    #ax.set_title(npzDataFile + ' eclipse %d\nFitted Min Of Total Light Curve = %.13f\nMin of Light Curve = %.13f\n"error" = Fit Min-Min = %.13f\nparameters = [%.10f, %.10f, %.10f, %.10f, %.10f]'%(EclipseNum,WDeclipseTime,jd[DataIndex_min],(WDeclipseTime-jd[DataIndex_min]),TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4]),size='small')
        ax.set_title(npzDataFile + ' eclipse %d\n'%(EclipseNum)+BMJDfile+ '\nFitted Min Of Total Light Curve = %.13f\nparameters = [%.10f, %.10f, %.10f, %.10f, %.10f]'%(WDeclipseTime,TFparams[0],TFparams[1],TFparams[2],TFparams[3],TFparams[4]),size='small')


        ax.plot(Xrange+jd[0],WD,'r', label = 'White Dwarf Light Curve Fitted',)
    else:
        FittedParams.append([np.nan])
        EclipseTimes.append(np.nan)
        ax.set_title(npzDataFile + ' eclipse %d\n'%(EclipseNum)+BMJDfile,size='small')
    ax.set_xlabel('MJD(TDB)',size='small')
    ax.set_ylabel('Photon Count Rate per Sec',size='small')


np.savez('/Scratch/dataProcessing/SDSS_J0926/AllDataBlue-fitresults.npz',EclipseTimes = EclipseTimes, FittedParams = FittedParams)
ax.legend(prop={'size':'small'})
plt.show()


