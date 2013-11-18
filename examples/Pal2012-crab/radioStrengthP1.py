from matplotlib import rcParams, rc
from spuriousRadioProbP1 import probsOfGRP

# common setup for matplotlib
params = {'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 14,
          'lines.linewidth':1.5,
          'text.fontsize': 14,
          'legend.fontsize': 16,
          'xtick.labelsize': 14,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 14}
# use of Sans Serif also in math mode
rc('text.latex', preamble='\usepackage{sfmath}')

rcParams.update(params)

import matplotlib.pyplot as plt

def cm2inch(cm):
    """Centimeters to inches"""
    return cm *0.393701

import numpy as np
import matplotlib
import scipy.stats
import tables

def nSigma(pvalue):
    return scipy.special.erfinv(pvalue)*np.sqrt(2.)
np.seterr(divide='ignore')
np.set_printoptions(threshold=np.nan)

radioCalibrationFactor = 162.
path = '/Scratch/dataProcessing/crabData2/'
nIdxToCheck = 81
nSigmaRadioCutoff = 3
nBins = 250
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_indices.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_skyW.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_spec5.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_specAll.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_KSswap.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
dataFilePath = path+'indPulseProfiles_{}sigma_P1_weak.h5'.format(nSigmaRadioCutoff)
dataFile = tables.openFile(dataFilePath,mode='r')
radioMax = dataFile.root.radioMax.read()
radioMean = dataFile.root.radioMean.read()
counts = dataFile.root.counts.read()
giantPulseNumbers = dataFile.root.giantPulseNumbers.read()
pulseNumberTable = dataFile.root.pulseNumberTable.read()
giantPulseNumberMask = dataFile.root.giantPulseNumberMask.read()
idxOffsets = dataFile.root.idxOffsets.read()
indProfiles = dataFile.root.indProfiles.read()
radioIndices = dataFile.root.radioIndices.read()

overlapPNs = np.load('overlapP1.npz')['overlap']
mainPulseMask = np.logical_not(np.in1d(giantPulseNumbers,overlapPNs))

radioMax = radioMax[mainPulseMask]
counts = counts[mainPulseMask]
giantPulseNumbers = giantPulseNumbers[mainPulseMask]
pulseNumberTable = pulseNumberTable[mainPulseMask]
giantPulseNumberMask = giantPulseNumberMask[mainPulseMask]
indProfiles = indProfiles[mainPulseMask]
radioIndices = radioIndices[mainPulseMask]

radioStrength = radioMax
print 'radioMax',np.min(radioMax),np.median(radioMax),np.max(radioMax)
print np.shape(radioMax,),np.shape(radioIndices),np.shape(counts)

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
print 'dimMask',np.sum(dimMask),np.shape(dimMask)

nRadioBins=14

probDict = probsOfGRP(nRadioBins) #prints out and plots values
print 'interp range',probDict['grpProbFuncRange']
radioBins = probDict['radioBins']
#probGRP = probDict['probGRP']
radioPeakMask = probDict['radioPeakMask']
grpProbFunc = probDict['grpProbFunc']
lowerFuncLimit,upperFuncLimit = probDict['grpProbFuncRange']
nPoints = 1000
filledStrengthBins = np.linspace(lowerFuncLimit,upperFuncLimit,nPoints)
probGRP = grpProbFunc(filledStrengthBins)
print 'radioPeakMask',np.sum(radioPeakMask),np.shape(radioPeakMask)

indProfilesMask = np.tile(giantPulseNumberMask,(np.shape(indProfiles)[2],1,1))
indProfilesMask = np.swapaxes(indProfilesMask,0,2)
indProfilesMask = np.swapaxes(indProfilesMask,0,1)
indProfilesMasked = np.ma.array(indProfiles,mask=indProfilesMask)

fig = plt.figure()
ax = fig.add_subplot(111)
handleMatshow = ax.matshow(counts[dimMask])
ax.set_aspect(1.0*len(idxOffsets)/len(radioMax))
fig.colorbar(handleMatshow)

nIdxOffsets = len(idxOffsets)
#goodPulseMask = np.logical_and(counts > 0,)
#sum over GRP index, to get number of nonzero pulses in each index
# this will be used to scale later
nPulsesPerIdx = np.array(np.sum(giantPulseNumberMask,axis=0),dtype=np.double).reshape((-1,1))
idx0 = np.searchsorted(idxOffsets,0)

cmap = matplotlib.cm.jet
histStart = 0.
histEnd = 1.
nBins=np.shape(indProfiles)[2]

_,phaseBinEdges = np.histogram(np.array([0]),range=(histStart,histEnd),bins=nBins)
grpProfile = np.ma.mean(indProfilesMasked.data[:,idx0],axis=0)
#surroundingProfiles = np.ma.mean(indProfilesMasked,axis=0)#average over iGRP
#avgSurroundingProfile = np.ma.mean(surroundingProfiles,axis=0)#average over iIdxOffset
#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(phaseBinEdges[0:-1],grpProfile,c='k')
#for iIdxOffset,idxOffset in enumerate(idxOffsets):
#    profile = surroundingProfiles[iIdxOffset]
#    idxFraction = (iIdxOffset+1.)/nIdxOffsets
#    color=cmap(idxFraction)
#    if idxOffset != 0:
#        ax.plot(phaseBinEdges[0:-1],profile,c=color)
    
#nBinsInPeak = 5
peakIdx = np.argmax(grpProfile)
#peakBins = range(peakIdx-nBinsInPeak//2,peakIdx+nBinsInPeak//2+nBinsInPeak%2)
peakBins = range(peakIdx-1,peakIdx+2)
giantPeakHeights = np.sum(indProfiles[:,idx0,peakBins],axis=1)
peakHeights = np.sum(indProfiles[:,:,peakBins],axis=2)
#index peakHeights[iGRP,iIdxOffset]

maskedPeakHeights = np.ma.array(peakHeights,mask=giantPulseNumberMask)
avgPeakHeights  = np.ma.mean(maskedPeakHeights,axis=1)#average over iIdxOffset i.e. average of surrounding pulses for each iGRP
giantRatio = giantPeakHeights/avgPeakHeights
opticalEnhancement = (giantPeakHeights-avgPeakHeights)/avgPeakHeights


radioStrengthCutoff = 0.155
radioCutoffMask = radioStrength >= radioStrengthCutoff
strongMask = np.logical_and(radioCutoffMask,dimMask)
print 'strongMask',np.sum(strongMask),np.shape(strongMask)
print 'radioCutoffMask',np.sum(radioCutoffMask),np.shape(radioCutoffMask)
finalMask = np.logical_and(strongMask,radioPeakMask)
print 'finalMask',np.sum(finalMask),np.shape(finalMask)
#finalMask = strongMask

fig = plt.figure()
ax = fig.add_subplot(111)
handleMatshow = ax.matshow(counts[finalMask])
ax.set_aspect(1.0*len(idxOffsets)/len(counts[finalMask]))
fig.colorbar(handleMatshow)

binProfile = np.mean(indProfiles[finalMask,idx0,:],axis=0)
surroundingProfiles = np.ma.mean(indProfilesMasked[finalMask,:],axis=0)
avgProfile = np.ma.mean(surroundingProfiles,axis=0)
minProfileIndex = np.argmin(avgProfile)
skyLevel = np.mean(avgProfile[minProfileIndex-3:minProfileIndex+3])#change to average of low points
avgProfileErrors = np.ma.std(surroundingProfiles,axis=0)/np.sqrt(nIdxOffsets)#std over iIdxOffset /sqrt(N) to get error in avgProfile
#add errors in quadrature
skySigma = np.sqrt(np.sum(avgProfileErrors[minProfileIndex-3:minProfileIndex+3]**2.))
print 'sky level',skyLevel,'+/-',skySigma
binProfile-=skyLevel
surroundingProfiles-=skyLevel
avgProfile-=skyLevel
indProfiles-=skyLevel

overallCoincidentProfile=binProfile
avgOverallProfile = avgProfile
stdProfile = np.ma.std(surroundingProfiles,axis=0)#std over iIdxOffset
stdProfile = np.sqrt(stdProfile**2+skySigma**2)#add sky error in quadrature

giantPeakHeight = np.sum(binProfile[peakBins])
peakHeight = np.sum(avgProfile[peakBins])
peakSigma = np.sqrt(np.sum(stdProfile[peakBins]**2))
enhancement = (giantPeakHeight-peakHeight)/peakHeight
overallEnhancement = enhancement
enhancementNSigma = (giantPeakHeight-peakHeight)/peakSigma
enhancementError = peakSigma/peakHeight
print 'peak enhancement of avg above',radioStrengthCutoff,':',enhancement,'+/-',enhancementError,'(',enhancementNSigma,' sigma)'
#probGRP = np.concatenate((probGRP,np.ones(nRadioBins-len(probGRP))))

#hist,radioBins = np.histogram(radioStrength,bins=radioBins)#bins=nRadioMeanBins)
radioBinned = np.digitize(radioStrength,bins=radioBins)
#fig = plt.figure()
#ax = fig.add_subplot(111)
enhancements = []
enhancementNSigmas = []
enhancementErrors = []
radioBinExpectedCenters = []
radioBinStd = []
for iBin,bin in enumerate(radioBins[0:-1]):
    binMask = radioBinned==(iBin+1)
    binMask = np.logical_and(binMask,radioPeakMask)
    binMask = np.logical_and(binMask,dimMask)
    print bin,np.sum(binMask)

    radioStrengthsInBin = radioStrength[binMask]
    radioBinExpectedCenters.append(np.mean(radioStrengthsInBin))
    radioBinStd.append(np.std(radioStrengthsInBin))
    
    binProfile = np.mean(indProfiles[binMask,idx0,:],axis=0)
    surroundingProfiles = np.ma.mean(indProfilesMasked[binMask,:],axis=0)
    avgProfile = np.ma.mean(surroundingProfiles,axis=0)
    stdProfile = np.ma.std(surroundingProfiles,axis=0)#std over iIdxOffset

    giantPeakHeight = np.sum(binProfile[peakBins])
    peakHeight = np.sum(avgProfile[peakBins])
    peakSigma = np.sqrt(np.sum(stdProfile[peakBins]**2))
    enhancement = (giantPeakHeight-peakHeight)/peakHeight
    enhancementNSigma = (giantPeakHeight-peakHeight)/peakSigma
    enhancementError = peakSigma/peakHeight
    enhancements.append(enhancement)
    enhancementNSigmas.append(enhancementNSigma)
    enhancementErrors.append(enhancementError)

    iBinFraction = (iBin+1.)/len(radioBins)
    colorFraction = iBinFraction
    color=cmap(colorFraction)
    #ax.plot(phaseBinEdges[0:-1],binProfile,c=color,label='%.3f-%.3f'%(radioBins[iBin],radioBins[iBin+1]))
    

enhancements = np.array(enhancements)
enhancementErrors = np.array(enhancementErrors)
radioBinExpectedCenters = np.array(radioBinExpectedCenters)
radioBinStd = np.array(radioBinStd)*radioCalibrationFactor
#Do the same for an all GRP above a threshold


fig = plt.figure()
ax = fig.add_subplot(111)
radioBinCenters = radioBins[0:-1]+np.diff(radioBins)/2.
#ax.errorbar(radioCalibrationFactor*radioBinCenters,enhancements,yerr=enhancementErrors,marker='.',color='g',label='Enhancement')
ax.errorbar(radioCalibrationFactor*radioBinCenters,100.*enhancements,yerr=100.*enhancementErrors,marker='.',color='k',linestyle='')
ax.plot(radioCalibrationFactor*np.array(radioBins),100.*np.append(enhancements,enhancements[-1]),label='Enhancement',drawstyle='steps-post',color='k')
ax.plot(radioCalibrationFactor*filledStrengthBins,100.*probGRP*overallEnhancement,color='.5',label='Enhancement depressed by\nFalse GRP Triggers',linestyle='--')
ax.set_xlabel('GRP Peak Flux Density (Jy)')
ax.set_xscale('log')
ax.set_ylabel('Optical Enhancement of GRP-Coincident Pulses (%)')
startTick = 18.
endTick = 28.
tickWidth = 1.
nTicks = (endTick-startTick)/tickWidth
print 'filledStrengthBins',filledStrengthBins[0],filledStrengthBins[-1]
print 'probGRP',100.*probGRP[0]*overallEnhancement
print startTick,endTick,tickWidth,nTicks
ticks = np.arange(startTick,endTick+tickWidth,tickWidth)
ax.set_xticks(ticks)
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%f'))
ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, pos: str(x)))
ax.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(2.))
ax.set_xlim([startTick,endTick])
ax.set_ylim((-1.5,7))
ax.legend(loc='best')

print 'range',radioBinCenters[0],radioBinCenters[-1]
smoothModel = np.array(grpProbFunc(radioBinCenters)*overallEnhancement)

chi2=np.sum((enhancements-smoothModel)**2/enhancementErrors**2)
dof=len(enhancements)-1
smoothPvalue=1-scipy.stats.chi2.cdf(chi2,dof)
print 'smooth chi2 dof pvalue',chi2,dof,smoothPvalue,nSigma(1-smoothPvalue),'sigmas'

plt.show()


