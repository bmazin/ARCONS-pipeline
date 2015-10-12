import numpy as np
import scipy.stats,scipy.special
import matplotlib.pyplot as plt
from util.popup import plotArray,PopUp
import scipy.ndimage
import os
from util.readDict import readDict

def plotPulseProfile(ax,phaseBinEdges,pulseProfile,profileErrors=None,plotDoublePulse=True,**kwargs):
    if plotDoublePulse:
        doublePhaseBinEdges = np.concatenate([phaseBinEdges,phaseBinEdges[1:]+1.])
        doubleSteppedPulseProfile = np.concatenate([pulseProfile,pulseProfile,[pulseProfile[-1]]])
        ax.plot(doublePhaseBinEdges,doubleSteppedPulseProfile,drawstyle='steps-post',**kwargs)
        if not (profileErrors is None):
            doublePulseProfile = np.concatenate([pulseProfile,pulseProfile])
            doubleProfileErrors = np.concatenate([profileErrors,profileErrors])
            doubleBinCenters = doublePhaseBinEdges[0:-1]+np.diff(doublePhaseBinEdges)/2.
            ax.errorbar(doubleBinCenters,doublePulseProfile,yerr=doubleProfileErrors,linestyle='',**kwargs)
    else:
        steppedPulseProfile = np.concatenate([pulseProfile,[pulseProfile[-1]]])
        ax.plot(phaseBinEdges,steppedPulseProfile,drawstyle='steps-post',**kwargs)
        if not (profileErrors is None):
            binCenters = phaseBinEdges[0:-1]+np.diff(phaseBinEdges)/2.
            ax.errorbar(binCenters,pulseProfile,yerr=profileErrors,linestyle='',**kwargs)

def nSigma(pvalue):
    return scipy.special.erfinv(pvalue)*np.sqrt(2.)

def weightedStd(values, **kwargs):
    """
    Return the weighted average and standard deviation.
    values, kwargs -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, **kwargs)
    variance = np.average((values-average)**2, **kwargs)  # Fast and numerically precise
    return {'avg':average,'std':np.sqrt(variance)}

def rebinHists(oldBinEdges,oldHist,nOldBinsInNewBin=2,axis=-1):
    firstAfterConvolve = (nOldBinsInNewBin//2) - 1 +(nOldBinsInNewBin%2)
    rebinnedEdges = oldBinEdges[::nOldBinsInNewBin]
    convolvedHists = scipy.ndimage.filters.convolve1d(oldHist,np.ones(nOldBinsInNewBin),mode='constant',axis=axis)
    rebinnedHists = convolvedHists[:,firstAfterConvolve::nOldBinsInNewBin]
    return {'newBinEdges':rebinnedEdges,'newHists':rebinnedHists}

def flatChiTest(phaseProfile,profileErrors,verbose=True):

    profileAvg = np.average(phaseProfile,weights=1./profileErrors**2)
    errorProfileAvg = 1./np.sqrt(np.sum(1./profileErrors**2))
    chi2=np.sum((phaseProfile-profileAvg)**2/profileErrors**2)
    dof = len(phaseProfile) - 1 #1 free parameter - the average
    reducedChi2 = 1.*chi2/dof
    pvalue=1-scipy.stats.chi2.cdf(chi2,dof)
    nSigmaSignificance = nSigma(1-pvalue)
    if verbose:
        print 'chi2:', chi2
        print 'reduced chi2:', reducedChi2
        print 'pvalue:',pvalue
        print 'significance:',nSigmaSignificance,'sigmas'
        print ''
    return {'profileAvg':profileAvg,'errorProfileAvg':errorProfileAvg,'chi2':chi2,'dof':dof,'reducedChi2':reducedChi2,'pvalue':pvalue,'nSigmaSignificance':nSigmaSignificance}

def combineProfiles(phaseBinEdges,phaseProfiles,nNewBins=10,label='',verbose=True,bPlotProfile=False):
    nOldBins = len(phaseBinEdges)-1
    rebinFactor = nOldBins//nNewBins
    rebinDict = rebinHists(phaseBinEdges,phaseProfiles,nOldBinsInNewBin=rebinFactor)
    phaseBinEdges = rebinDict['newBinEdges']
    phaseProfiles = rebinDict['newHists']

    phaseProfile = np.sum(phaseProfiles,axis=0)

    allErrors = np.std(phaseProfiles)
    #profileErrors = np.std(phaseProfiles,axis=0)
    profileErrors = np.sqrt(phaseProfile)


    if verbose:
        print label,'analysis'
        print 'weighted average(%):',profileAvg,'+/-',errorProfileAvg

    if bPlotProfile:
        fig,ax = plt.subplots()
        plotPulseProfile(ax,phaseBinEdges,phaseProfile,profileErrors,color='b',plotDoublePulse=False)
        ax.set_xlabel('phase')
        ax.set_ylabel('counts')
        ax.set_title('J0337 Pulse Profile - {}'.format(label))

    return {'phaseBinEdges':phaseBinEdges,'phaseProfile':phaseProfile,'profileErrors':profileErrors}

if __name__=='__main__':


    paramFile = 'j0337.dict'
    params = readDict()
    params.read_from_file(paramFile)
    obsSequences = [params[key] for key in sorted(params.keys()) if key.startswith('obsSequence')]
    sunsetDates = [params[key] for key in sorted(params.keys()) if key.startswith('sunsetDate')]

    tstampsToUse = np.concatenate(obsSequences)

    nPhaseBins = 150
    wvlStart = 4000 #angstrom
    wvlEnd = 5500 #angstrom
    #dataPath = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_{}arcsecAperture.npz'.format(nPhaseBins,wvlStart,wvlEnd,apertureRadius)
    dataPath = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_optimalAperture.npz'.format(nPhaseBins,wvlStart,wvlEnd)
    print os.path.basename(dataPath)

    dataDict = np.load(dataPath)
    phaseProfiles = np.array(dataDict['phaseProfiles'],dtype=np.double)
    phaseBinEdges = dataDict['phaseBinEdges']
    tstamps = dataDict['obsTimestamps']
    apertureRadiusList = dataDict['apertureRadiusList']
    psfFits = dataDict['psfFits']
    print len(tstamps),'files in phaseProfiles'

    tstampMask = np.array([tstamp in tstampsToUse for tstamp in tstamps])
    print np.sum(tstampMask), 'files used'
    tstamps = tstamps[tstampMask]
    phaseProfiles = phaseProfiles[tstampMask]
    apertureRadiusList = apertureRadiusList[tstampMask]
    psfFits = psfFits[tstampMask]

    combineDict = combineProfiles(phaseBinEdges,phaseProfiles,label='total',verbose=False,bPlotProfile=False)
    #return {'profileAvg':profileAvg,'errorProfileAvg':errorProfileAvg,'chi2':chi2,'dof':dof,'reducedChi2':reducedChi2,'pvalue':pvalue,'nSigmaSignificance':nSigmaSignificance}

    chiDict = flatChiTest(combineDict['phaseProfile'],combineDict['profileErrors'],verbose=True)
    print 'avg',chiDict['profileAvg']
    phaseBinEdges = combineDict['phaseBinEdges']

    nSigmaThreshold = 5.

    #sharp pulsar

    normSharpPulseProfile = np.zeros_like(combineDict['phaseProfile'])
    fakePulseIdx = 2
    normSharpPulseProfile[fakePulseIdx] = 1.

    nSigmaSig = 0.
    fakePulseAmplitude = 100.
    while  nSigmaSig < nSigmaThreshold:

        fakedProfile = combineDict['phaseProfile'] + fakePulseAmplitude*normSharpPulseProfile
        fakedProfileErrors = np.sqrt(fakedProfile)

        fakedChiDict = flatChiTest(fakedProfile,fakedProfileErrors,verbose=False)
        nSigmaSig = fakedChiDict['nSigmaSignificance']
        if nSigmaSig < nSigmaThreshold:
            fakePulseAmplitude += 100.
    sharpAmp = fakePulseAmplitude

    print 'sharp peaked profile'
    print 'amplitude',sharpAmp
    sharpCps = np.mean(sharpAmp*normSharpPulseProfile)
    print 'cps',sharpCps
    fluxFraction = 1.*sharpAmp/chiDict['profileAvg']
    print 'fraction',fluxFraction
    print 'mag diff',-2.5*np.log10(fluxFraction)
    print 'phase averaged',-2.5*np.log10(1.*sharpCps/chiDict['profileAvg'])

    
    #double peak - B0656
    normDoublePulseProfile = np.array([
    -0.04261796042617849,
    0.49923896499239007,
    2.532724505327245,
    1.442922374429224,
    1.1567732115677325,
    0.45053272450532766,
    0.3774733637747345,
    1.1993911719939119,
    2.28310502283105,
    0.14003044140030507])

    normDoublePulseProfile -= np.min(normDoublePulseProfile)

    nSigmaSig = 0.
    fakePulseAmplitude = 100.
    while  nSigmaSig < nSigmaThreshold:
        fakedProfile = combineDict['phaseProfile'] + fakePulseAmplitude* normDoublePulseProfile
        fakedProfileErrors = np.sqrt(fakedProfile)

        fakedChiDict = flatChiTest(fakedProfile,fakedProfileErrors,verbose=False)
        nSigmaSig = fakedChiDict['nSigmaSignificance']
        if nSigmaSig < nSigmaThreshold:
            fakePulseAmplitude += 100.

    doubleAmp = fakePulseAmplitude
    flux = np.mean(doubleAmp*normDoublePulseProfile)
    print ''
    print 'double peaked profile'
    print 'amplitude',doubleAmp, 'flux',flux
    fluxFraction = 1.*flux/chiDict['profileAvg']
    print 'fraction',fluxFraction
    print 'mag diff',-2.5*np.log10(fluxFraction)

    nSigmaSig = 0.
    fakePulseAmplitude = 100.
    normSinePulseProfile = np.zeros_like(combineDict['phaseProfile'])
    phaseBinCenters = phaseBinEdges[0:-1] + np.diff(phaseBinEdges)/2.
    normSinePulseProfile = np.sin(2.*np.pi*phaseBinCenters)
    normSinePulseProfile -= np.min(normSinePulseProfile)
    nomrSinePulseProfile = 1.*normSinePulseProfile / np.sum(normSinePulseProfile)

    while  nSigmaSig < nSigmaThreshold:
        fakedProfile = combineDict['phaseProfile'] + fakePulseAmplitude* normSinePulseProfile
        fakedProfileErrors = np.sqrt(fakedProfile)

        fakedChiDict = flatChiTest(fakedProfile,fakedProfileErrors,verbose=False)
        nSigmaSig = fakedChiDict['nSigmaSignificance']
        if nSigmaSig < nSigmaThreshold:
            fakePulseAmplitude += 100.

    sineAmp = fakePulseAmplitude
    flux = np.mean(sineAmp*normSinePulseProfile)
    totalSineProfile = np.array(fakedProfile)
    print ''
    print 'sine profile'
    print 'amplitude',sineAmp, 'flux',flux, 'photons',np.sum(sineAmp*normSinePulseProfile)
    fluxFraction = 1.*flux/chiDict['profileAvg']
    print 'fraction',fluxFraction
    print 'mag diff',-2.5*np.log10(fluxFraction)

    fig,axs = plt.subplots(4,1)

    plotPulseProfile(axs[0],phaseBinEdges,combineDict['phaseProfile'],combineDict['profileErrors'],color='black',plotDoublePulse=False)
    plotPulseProfile(axs[1],phaseBinEdges,sharpAmp*normSharpPulseProfile,np.sqrt(sharpAmp*normSharpPulseProfile),color='b',plotDoublePulse=False)
    plotPulseProfile(axs[2],phaseBinEdges,doubleAmp*normDoublePulseProfile,np.sqrt(doubleAmp*normDoublePulseProfile),color='r',plotDoublePulse=False)

    plotPulseProfile(axs[3],phaseBinEdges,combineDict['phaseProfile'],combineDict['profileErrors'],color='black',plotDoublePulse=False)
    plotPulseProfile(axs[3],phaseBinEdges,combineDict['phaseProfile']+sharpAmp*normSharpPulseProfile,np.sqrt(combineDict['phaseProfile']+sharpAmp*normSharpPulseProfile),color='b',plotDoublePulse=False)
    plotPulseProfile(axs[3],phaseBinEdges,combineDict['phaseProfile']+doubleAmp*normDoublePulseProfile,np.sqrt(combineDict['phaseProfile']+doubleAmp*normDoublePulseProfile),color='r',plotDoublePulse=False)

    fig,ax = plt.subplots(1,1)
    plotPulseProfile(ax,phaseBinEdges,combineDict['phaseProfile'],combineDict['profileErrors'],color='black',plotDoublePulse=False)
    plotPulseProfile(ax,phaseBinEdges,totalSineProfile,np.sqrt(totalSineProfile),color='r',plotDoublePulse=False)
    

    plt.show()

    

