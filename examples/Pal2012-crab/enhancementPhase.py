from matplotlib import rcParams, rc
from spuriousRadioProbRange import probsOfGRP
from util import mpfit
from util.fitFunctions import gaussian
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import scipy.stats
import tables
import scipy.special


def fitGauss(xdata,ydata,yerr,flatLine=False):

    nBins=100
    amplitude = .5*np.max(ydata)
    x_offset = xdata[np.argmax(ydata)]
    sigma = (np.max(xdata)-np.min(xdata))/10.
    y_offset = 3.
    fixed = [False]*4
    if flatLine == True:
        amplitude = 0
        fixed[0:3] = [True]*3

    params=[sigma, x_offset, amplitude, y_offset]  # First guess at fit params
    errs = yerr
    errs[np.where(errs == 0.)] = 1.
    quiet = True

    parinfo = [ {'n':0,'value':params[0],'limits':[.0001, .1], 'limited':[True,True],'fixed':fixed[0],'parname':"Sigma",'error':0},
       {'n':1,'value':params[1],'limits':[x_offset-sigma*3, x_offset+sigma*3],'limited':[True,True],'fixed':fixed[1],'parname':"x offset",'error':0},
       {'n':2,'value':params[2],'limits':[.2*amplitude, 3.*amplitude],'limited':[True,True],'fixed':fixed[2],'parname':"Amplitude",'error':0},
       {'n':3,'value':params[3],'limited':[False,False],'fixed':fixed[3],'parname':"y_offset",'error':0}]

    fa = {'x':xdata,'y':ydata,'err':yerr}

    m = mpfit.mpfit(gaussian, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
    if m.status <= 0:
        print m.status, m.errmsg

    mpp = m.params                                #The fit params
    mpperr = m.perror

    for k,p in enumerate(mpp):
        parinfo[k]['value'] = p
        parinfo[k]['error'] = mpperr[k]
        #print parinfo[k]['parname'],p," +/- ",mpperr[j]
        if k==0: sigma = p
        if k==1: x_offset = p
        if k==2: amplitude = p
        if k==3: y_offset = p

    fineXdata = np.linspace(np.min(xdata),np.max(xdata),100.)
    gaussfit = y_offset + amplitude * np.exp( - (( xdata - x_offset)**2) / ( 2. * (sigma**2)))
    fineGaussFit = y_offset + amplitude * np.exp( - (( fineXdata - x_offset)**2) / ( 2. * (sigma**2)))

    resolution = np.abs(x_offset/(2.355*sigma))
    return {'gaussfit':gaussfit,'resolution':resolution,'sigma':sigma,'x_offset':x_offset,'amplitude':amplitude,'y_offset':y_offset,'fineXdata':fineXdata,'fineGaussFit':fineGaussFit,'parinfo':parinfo}



# common setup for matplotlib
params = {'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 14,
          'lines.linewidth': 1.5,
          'text.fontsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 14,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 14}
# use of Sans Serif also in math mode
rc('text.latex', preamble='\usepackage{sfmath}')

rcParams.update(params)

phaseShift = 1.-0.677001953125#found with findOpticalPeak.py

def align_yaxis(ax1, v1, ax2, v2):
    """
    adjust ax2 ylimit so that v2 in ax2 is aligned to v1 in ax1
    Taken from http://stackoverflow.com/questions/10481990/matplotlib-axis-with-two-scales-shared-origin
    """
    _, y1 = ax1.transData.transform((0, v1))
    _, y2 = ax2.transData.transform((0, v2))
    inv = ax2.transData.inverted()
    _, dy = inv.transform((0, 0)) - inv.transform((0, y1-y2))
    miny, maxy = ax2.get_ylim()
    ax2.set_ylim(miny+dy, maxy+dy)

def indexToPhase(indices,radioIndexOffset=0.5):
    radioArrivalPhases = (indices+radioIndexOffset)/2048.+phaseShift
    return radioArrivalPhases

def nSigma(pvalue):
    return scipy.special.erfinv(pvalue)*np.sqrt(2.)

np.seterr(divide='ignore')
np.set_printoptions(threshold=np.nan)

path = '/Scratch/dataProcessing/crabData2/'
nIdxToCheck = 81
nSigmaRadioCutoff = 3
nBins = 250
bUseFineIndexBins = False
bInterpulses = False

#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_swap.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
dataFilePath = path+'indPulseProfiles_{}sigma_P2_KS.h5'.format(nSigmaRadioCutoff)
dataFile = tables.openFile(dataFilePath,mode='r')
radioMax = dataFile.root.radioMax.read()
counts = dataFile.root.counts.read()#-dataFile.root.skyCounts.read()
giantPulseNumbers = dataFile.root.giantPulseNumbers.read()
pulseNumberTable = dataFile.root.pulseNumberTable.read()
giantPulseNumberMask = dataFile.root.giantPulseNumberMask.read()
idxOffsets = dataFile.root.idxOffsets.read()
indProfiles = dataFile.root.indProfiles.read()
radioIndices = dataFile.root.radioIndices.read()

radioIndexBins=np.array([1369,1371,1373,1375,1378,1381,1385,1389,1395])-.5
radioIndexBinsFine = np.arange(1369,1396)-.5

if bUseFineIndexBins == True:
    radioIndexBins = radioIndexBinsFine


startRadioIndex = radioIndexBins[0]
endRadioIndex = radioIndexBins[-1]

probDict = probsOfGRP(startPeakIndex=startRadioIndex,endPeakIndex=endRadioIndex)
probPhaseBins = probDict['radioPhaseBins']
probPeakDist = probDict['peakDist']


#a mask for less good data, during bright or dim times
dimMask = np.ones(len(counts))
idx0 = np.searchsorted(idxOffsets,0)
dimMask[counts[:,idx0]==0]=0
lineCounts = np.mean(counts,axis=1)

meanLineCounts = np.mean(lineCounts[lineCounts!=0])
stdLineCounts = np.std(lineCounts[lineCounts!=0])
stdPercentCutoff=0.
upperCutoff = scipy.stats.scoreatpercentile(lineCounts,100.-stdPercentCutoff)
lowerCutoff = scipy.stats.scoreatpercentile(lineCounts,stdPercentCutoff)
dimMask[lineCounts>upperCutoff] = 0
dimMask[lineCounts<lowerCutoff] = 0
dimMask = (dimMask==1)

radioStrength = radioMax



indProfilesMask = np.tile(giantPulseNumberMask,(np.shape(indProfiles)[2],1,1))
indProfilesMask = np.swapaxes(indProfilesMask,0,2)
indProfilesMask = np.swapaxes(indProfilesMask,0,1)
indProfilesMasked = np.ma.array(indProfiles,mask=indProfilesMask)

nIdxOffsets = len(idxOffsets)
#sum over GRP index, to get number of nonzero pulses in each index
# this will be used to scale later
nPulsesPerIdx = np.array(np.sum(giantPulseNumberMask,axis=0),dtype=np.double).reshape((-1,1))

cmap = matplotlib.cm.jet
histStart = 0.
histEnd = 1.
nBins=np.shape(indProfiles)[2]

_,phaseBinEdges = np.histogram(np.array([0]),range=(histStart,histEnd),bins=nBins)
phaseBinEdges+=phaseShift
phaseBinCenters = phaseBinEdges[0:-1]+np.diff(phaseBinEdges)/2.
grpProfile = np.ma.mean(indProfilesMasked.data[:,idx0],axis=0)
    
peakIdx = np.argmax(grpProfile)
peakBins = range(peakIdx-1,peakIdx+2)
print 'opticalPeakPhaseBins',peakBins




nRadioBins=15

radioStrengthCutoff = .175#0.175
radioCutoffMask = radioStrength >= radioStrengthCutoff
strongMask = np.logical_and(radioCutoffMask,dimMask)
#finalMask = np.logical_and(strongMask,radioPeakMask)
radioPhaseMask = np.logical_and(radioIndices >= 1369,radioIndices <= 1394)
finalMask = np.logical_and(strongMask,radioPhaseMask)
print 'GRP above',radioStrengthCutoff,':',np.sum(finalMask),'and in phase range'

#counts color plot
fig = plt.figure()
ax = fig.add_subplot(111)
handleMatshow = ax.matshow(counts[finalMask])
ax.set_aspect(1.0*np.shape(counts[finalMask])[1]/np.shape(counts[finalMask])[0])
fig.colorbar(handleMatshow)

overallCoincidentProfile = np.mean(indProfiles[finalMask,idx0,:],axis=0)
surroundingProfiles = np.ma.mean(indProfilesMasked[finalMask,:],axis=0)
avgProfile = np.ma.mean(surroundingProfiles,axis=0)
minProfileIndex = np.argmin(avgProfile)

#for the sky level take an average over 5 points at the lowest part of the period
skyLevel = np.mean(avgProfile[minProfileIndex-3:minProfileIndex+3])
avgProfileErrors = np.ma.std(surroundingProfiles,axis=0)/np.sqrt(nIdxOffsets)#std over iIdxOffset /sqrt(N) to get error in avgProfile
#add errors in quadrature
skySigma = np.sqrt(np.sum(avgProfileErrors[minProfileIndex-3:minProfileIndex+3]**2.))
#should check error in sky level at some point
print 'sky level',skyLevel,'+/-',skySigma
overallCoincidentProfile-=skyLevel
surroundingProfiles-=skyLevel
avgProfile-=skyLevel
indProfiles-=skyLevel
avgOverallProfile = avgProfile
stdProfile = np.ma.std(surroundingProfiles,axis=0)#std over iIdxOffset
stdProfile = np.sqrt(stdProfile**2+skySigma**2)
avgStdProfile = stdProfile/np.sqrt(nIdxOffsets-1)


giantPeakHeight = np.sum(overallCoincidentProfile[peakBins])
peakHeight = np.sum(avgProfile[peakBins])
peakSigma = np.sqrt(np.sum(stdProfile[peakBins]**2))
overallEnhancement = (giantPeakHeight-peakHeight)/peakHeight
enhancementNSigma = (giantPeakHeight-peakHeight)/peakSigma
enhancementError = peakSigma/peakHeight
overallEnhancementError = enhancementError
print 'peak enhancement of avg above',radioStrengthCutoff,':',overallEnhancement,'+/-',enhancementError,'(',enhancementNSigma,' sigma)'

overallPeakHeight = np.array(peakHeight)

allProfiles = np.array(surroundingProfiles.data)
allProfiles[idx0]=overallCoincidentProfile#add back in since it was masked and zeroed earlier
allPeakHeights = np.sum(allProfiles[:,peakBins],axis=1)
peakPercentDifferenceByIdxOffset = (allPeakHeights-peakHeight)/peakHeight 
nSigmaByIdxOffset = (allPeakHeights-peakHeight)/peakSigma

#significance figure
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(idxOffsets,np.abs(nSigmaByIdxOffset),'k')
ax.set_ylabel('Standard Deviations of Peak Height from Average Peak')
ax.set_xlabel('Pulse Offset Relative to GRP (number of periods)')
ax.set_ylim((0,8))

np.savez('sigP2.npz',idxOffsets=idxOffsets,nSigmaByIdxOffset=nSigmaByIdxOffset)

giantPeakHeights = np.sum(indProfiles[:,idx0,peakBins][finalMask],axis=1)
peakHeights = np.sum(indProfiles[:,:,peakBins][finalMask],axis=2)
#index peakHeights[iGRP,iIdxOffset]

maskedPeakHeights = np.ma.array(peakHeights,mask=giantPulseNumberMask[finalMask])
avgPeakHeights  = np.ma.mean(maskedPeakHeights,axis=1)#average over iIdxOffset i.e. average of surrounding pulses for each iGRP
opticalEnhancementGRP = (giantPeakHeights-avgPeakHeights)/avgPeakHeights
opticalEnhancement = (avgPeakHeights-overallPeakHeight)/overallPeakHeight


radioProfile = np.loadtxt(path+'radio/RadioProfile_LyneDM_TZRCorrect_withGUPPIdelay.txt',skiprows=1,usecols=[3])
nRadioPhaseBins = len(radioProfile)
radioProfilePhaseBins = (1.*np.arange(nRadioPhaseBins)+.5)/nRadioPhaseBins
radioProfilePhaseBins+=phaseShift

fig = plt.figure()
ax = fig.add_subplot(111)
ax2 = ax.twinx()
pltHandle2 = ax2.plot(radioProfilePhaseBins,radioProfile,c=(.4,.5,.8),label='Radio Pulse')
pltHandle0 = ax.errorbar(phaseBinCenters,overallCoincidentProfile,yerr=stdProfile,c='k',label='Optical GRP-accompanied Pulse')
pltHandle1 = ax.plot(phaseBinCenters,avgProfile,c='r',label='Optical non-GRP-accompanied Pulse')
pltHandles = [pltHandle0,pltHandle1[0],pltHandle2[0]]
pltLabels = [pltHandle.get_label() for pltHandle in pltHandles]

box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.85])
ax2.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 0.85])
ax.set_ylim((0.055,.081))
ax2.set_ylim((.11,.155))
ax.set_xlim((0.97,1.005))
locator = matplotlib.ticker.MultipleLocator(.01)
ax2.yaxis.set_major_locator(locator)
ax.legend(pltHandles,pltLabels,loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=2)
ax.set_ylabel('Optical Counts per Period per Pixel')
ax.set_xlabel('Phase')
ax2.set_ylabel('Normalized Radio Intensity')

#enhanced profile figure
#fig = plt.figure(figsize=(1.8,2))
ax = fig.add_subplot(2,2,1)
#ax = fig.add_axes([0.,.6,.4,.4])

doublePhaseBins = np.concatenate([phaseBinCenters-1,phaseBinCenters,1+phaseBinCenters])
doubleOverallCoincidentProfile = np.concatenate([overallCoincidentProfile,overallCoincidentProfile,overallCoincidentProfile])
doubleStdProfile = np.concatenate([stdProfile,stdProfile,stdProfile])
doubleAvgProfile = np.concatenate([avgProfile,avgProfile,avgProfile])
doubleRadioProfilePhaseBins = np.concatenate([radioProfilePhaseBins-1,radioProfilePhaseBins,1+radioProfilePhaseBins])
doubleRadioProfile = np.concatenate([radioProfile,radioProfile,radioProfile])
ax2 = ax.twinx()
pltHandle2 = ax2.plot(doubleRadioProfilePhaseBins,doubleRadioProfile,c=(.4,.5,.8),label='Radio Pulse')
pltHandle0 = ax.plot(doublePhaseBins,doubleOverallCoincidentProfile,c='k',label='Optical GRP-accompanied Pulse')
pltHandle1 = ax.plot(doublePhaseBins,doubleAvgProfile,c='r',label='Optical non-GRP-accompanied Pulse')
pltHandles = [pltHandle0[0],pltHandle1[0],pltHandle2[0]]
pltLabels = [pltHandle.get_label() for pltHandle in pltHandles]

#rect = plt.Rectangle((.970,0.055),1.005-.970,.081-0.055,edgecolor='green',fill=True,linewidth=2.)
#ax.add_patch(rect)

ax.yaxis.set_visible(False)
ax2.yaxis.set_visible(False)

ax.set_ylim((-.005,.081))
ax2.set_ylim((-.01,.155))
ax.set_xlim((0.01,1.99))
ax.set_xlabel('Phase')
#ax.xaxis.label.set_size(14)
#ax.tick_params(axis='both', which='major', labelsize=12)
#ax2.tick_params(axis='both', which='major', labelsize=12)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.3,
                 box.width*.8, box.height * 0.7])
ax2.set_position([box.x0, box.y0 + box.height * 0.3,
                 box.width*.8, box.height * 0.7])



radioProfile = radioProfile*np.max(overallCoincidentProfile)/np.max(radioProfile)


#Now plot optical enhancement vs radio arrival time

radioPhases = indexToPhase(radioIndices)
radioPhaseBins = indexToPhase(radioIndexBins)

radioPhaseBinCenters = radioPhaseBins[0:-1]+np.diff(radioPhaseBins)/2.

print 'radioIndexBins',radioIndexBins,np.diff(radioIndexBins)
print 'radioPhaseBinEdges',radioPhaseBins,np.diff(radioPhaseBins)

radioBinned = np.digitize(radioPhases,bins=radioPhaseBins)
enhancements = []
enhancementNSigmas = []
enhancementNSigmasOverOverall = []
enhancementErrors = []
globalEnhancements = []
globalEnhancementErrors = []
#non-GRP pulse enhancements
profiles = []
radioStrengthCutoff = .175#0.175
radioCutoffMask = radioStrength >= radioStrengthCutoff
strongMask = np.logical_and(radioCutoffMask,dimMask)
strongMask = np.logical_and(strongMask,radioPhaseMask)

for iBin,bin in enumerate(radioPhaseBins[0:-1]):
    binMask = np.logical_and(radioBinned==(iBin+1),strongMask)
    binProfile = np.mean(indProfiles[binMask,idx0,:],axis=0)
    profiles.append(binProfile)
    phasesInBin = radioPhases[binMask]

    surroundingProfiles = np.ma.mean(indProfilesMasked[binMask,:],axis=0)

    avgProfile = np.ma.mean(surroundingProfiles,axis=0)
    stdProfile = np.ma.std(surroundingProfiles,axis=0)#std over iIdxOffset
    nSurroundingProfiles=np.sum(np.logical_not(surroundingProfiles.mask),axis=0)
    errorAvgProfile = np.divide(stdProfile,np.sqrt(nSurroundingProfiles))

    giantPeakHeight = np.sum(binProfile[peakBins])
    peakHeight = np.sum(avgProfile[peakBins])
    peakSigma = np.sqrt(np.sum(stdProfile[peakBins]**2))
    enhancement = (giantPeakHeight-peakHeight)/peakHeight
    enhancementNSigma = (giantPeakHeight-peakHeight)/peakSigma
    enhancementError = peakSigma/peakHeight
    enhancements.append(enhancement)
    enhancementNSigmas.append(enhancementNSigma)
    enhancementErrors.append(enhancementError)
    nSigmaOverOverall=(enhancement-overallEnhancement)/enhancementError
    enhancementNSigmasOverOverall.append(nSigmaOverOverall)
    #print '{:.3}+/-{:.3}({:.3},{:.3})'.format(enhancement,enhancementError,enhancementNSigma,(enhancement-overallEnhancement)/enhancementError)
    print '{}\t{:.5}\t{}\t{:.5}\t{:.5}\t{:.5}'.format(radioIndexBins[iBin],bin,np.sum(binMask),enhancement,enhancementError,nSigmaOverOverall)

    globalEnhancement = (giantPeakHeight-overallPeakHeight)/overallPeakHeight
    globalEnhancementError = peakSigma/overallPeakHeight
    globalEnhancements.append(globalEnhancement)
    globalEnhancementErrors.append(globalEnhancementError)

    nonGRPEnhancement = (peakHeight-overallPeakHeight)/overallPeakHeight
    nonGRPPeakSigma = np.sqrt(np.sum(errorAvgProfile[peakBins])**2)
    nonGRPEnhancementNSigma = (peakHeight-overallPeakHeight)/nonGRPPeakSigma
    nonGRPEnhancementError = nonGRPPeakSigma/overallPeakHeight
    #print 'nonGRP {:.3}+/-{:.3}({:.3})'.format(nonGRPEnhancement,nonGRPEnhancementError,nonGRPEnhancementNSigma)

    nextBin = radioPhaseBins[iBin+1]
    #ax.plot(phaseBinEdges[0:-1],binProfile-avgProfile,c=color,label='{:.3}-{:.3}'.format(bin,nextBin))
    #ax2.errorbar(phaseBinEdges[0:-1],binProfile,yerr=stdProfile,c=color,label='{:.3}-{:.3}'.format(bin,nextBin))
    #ax3.errorbar(phaseBinEdges[0:-1],avgProfile,yerr=errorAvgProfile,c=color,label='{:.3}-{:.3}'.format(bin,nextBin))
enhancements = np.array(enhancements)
enhancementErrors = np.array(enhancementErrors)
    
percentEnhancements = 100.*enhancements
percentEnhancementErrors = 100.*enhancementErrors

fig = plt.figure(figsize=(8.,6.))
#ax = fig.add_subplot(211)
ax = fig.add_axes([.15,.6,.8,.3])
#ax.step(radioIndexBins[0:-1],noiseDist,'g',label='noise detections')
ax.plot(probPhaseBins,np.append(probPeakDist,probPeakDist[-1]),'k',drawstyle='steps-post',label='GRP+noise detections')
ax.set_ylabel('Number of\nGRPs detected')
ax.xaxis.set_visible(False)
ax.xaxis.set_ticks([])
#ax.step(radioIndexBins[0:-1],peakDist,'k',label='GRP+noise detections')

#fig = plt.figure()
#ax = fig.add_subplot(212)
ax2 = fig.add_axes([.15,.1,.8,.5])


#ax.errorbar(radioPhaseBinCenters,100.*enhancements,yerr=100.*enhancementErrors,marker='.',color='k',label='enhancement relative to surrounding nonGRP',linestyle='.')

ax2.errorbar(radioPhaseBinCenters,percentEnhancements,yerr=percentEnhancementErrors,linestyle='.',color='k')
ax2.plot(radioPhaseBins,np.append(percentEnhancements,percentEnhancements[-1]),'k',drawstyle='steps-post',label='enhancement relative to surrounding nonGRP')

opticalPeakPhase = 0.993998046875
ax2.axhline(0.,linewidth=1.,c='k')
ax2.axvline(opticalPeakPhase,c='gray',linestyle='--')

ax2.set_xlabel('GRP Arrival Phase')
ax2.set_ylabel('Optical Enhancement of\nGRP-Accompanied Pulses (%)')
ax2.set_ylim((-1.,16))
ax2.set_xlim((.991,1.005))
ax.set_xlim((.991,1.005))
fig.text(.175,.85,'(a)',size=16)
ax2.yaxis.get_major_ticks()[-1].label1.set_visible(False)
fig.text(.175,.55,'(b)',size=16)
#ax2.legend(loc='lower left')

fig = plt.figure()
ax = fig.add_subplot(111)

#ax.errorbar(radioPhaseBinCenters,100.*enhancements,yerr=100.*enhancementErrors,marker='.',color='k',label='enhancement relative to surrounding nonGRP',linestyle='.')

ax.errorbar(radioPhaseBinCenters,percentEnhancements,yerr=percentEnhancementErrors,linestyle='.',color='k')
ax.plot(radioPhaseBins,np.append(percentEnhancements,percentEnhancements[-1]),'k',drawstyle='steps-post',label='enhancement relative to surrounding nonGRP')
radioPhaseBinWidths = np.diff(radioPhaseBins)



meanPercentEnhancement = np.average(percentEnhancements,weights = 1/percentEnhancementErrors**2)
errorMeanPercentEnhancement = 1./np.sqrt(np.sum(1/percentEnhancementErrors**2))
print 'weighted average enhancement (%):',meanPercentEnhancement,'+/-',errorMeanPercentEnhancement
chi2=np.sum((percentEnhancements-meanPercentEnhancement)**2/percentEnhancementErrors**2)
dof=len(percentEnhancements)-1 #free parameter: meanEnhancement
pvalue=1-scipy.stats.chi2.cdf(chi2,dof)
print 'flat line: chi2 dof pvalue significance',chi2,dof,pvalue,nSigma(1-pvalue),'sigmas'

gaussDict = fitGauss(xdata=radioPhaseBinCenters,ydata=percentEnhancements,yerr=percentEnhancementErrors)
fit = gaussDict['gaussfit']
ax.plot(radioPhaseBinCenters,fit)
ax.plot(gaussDict['fineXdata'],gaussDict['fineGaussFit'])


chi2Fit=np.sum((percentEnhancements-fit)**2/percentEnhancementErrors**2)
dofFit=len(percentEnhancements)-4 #free parameters:y_offset,x_offset,amplitude,sigma
pvalueFit=1-scipy.stats.chi2.cdf(chi2Fit,dofFit)
print 'gaussian: chi2 dof pvalue significance',chi2Fit,dofFit,pvalueFit,nSigma(1-pvalueFit),'sigmas'
print gaussDict['parinfo']


flatLineDict = fitGauss(xdata=radioPhaseBinCenters,ydata=percentEnhancements,yerr=percentEnhancementErrors,flatLine=True)
fit = flatLineDict['gaussfit']
ax.plot(radioPhaseBinCenters,fit)
ax.plot(flatLineDict['fineXdata'],flatLineDict['fineGaussFit'])


chi2=np.sum((percentEnhancements-fit)**2/percentEnhancementErrors**2)
dof=len(percentEnhancements)-1 #free parameters:y_offset
pvalue=1-scipy.stats.chi2.cdf(chi2,dof)
print 'flatLine: chi2 dof pvalue significance',chi2,dof,pvalue,nSigma(1-pvalue),'sigmas'
print flatLineDict['parinfo']

chi2Diff = chi2-chi2Fit
dofDiff = dof-dofFit
pvalueDiff=1-scipy.stats.chi2.cdf(chi2Diff,dofDiff)
print 'diff: chi2 dof pvalue significance',chi2Diff,dofDiff,pvalueDiff,nSigma(1-pvalueDiff),'sigmas'





plt.show()


