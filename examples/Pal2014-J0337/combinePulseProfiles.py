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


def combineProfiles(phaseBinEdges,phaseProfiles,nNewBins=10,label='',ax=None):
    nOldBins = len(phaseBinEdges)-1
    rebinFactor = nOldBins//nNewBins
    rebinDict = rebinHists(phaseBinEdges,phaseProfiles,nOldBinsInNewBin=rebinFactor)
    phaseBinEdges = rebinDict['newBinEdges']
    phaseProfiles = rebinDict['newHists']

    phaseProfile = np.sum(phaseProfiles,axis=0)

    allErrors = np.std(phaseProfiles)
    #profileErrors = np.std(phaseProfiles,axis=0)
    profileErrors = np.sqrt(phaseProfile)

    profileAvg = np.average(phaseProfile,weights=1./profileErrors**2)
    errorProfileAvg = 1./np.sqrt(np.sum(1./profileErrors**2))

    print label,'analysis'
    print 'weighted average(%):',profileAvg,'+/-',errorProfileAvg
    chi2=np.sum((phaseProfile-profileAvg)**2/profileErrors**2)
    dof = len(phaseProfile) - 1 #1 free parameter - the average
    reducedChi2 = 1.*chi2/dof
    pvalue=1-scipy.stats.chi2.cdf(chi2,dof)
    print 'compare {} to flat line'.format(label)
    print 'chi2:', chi2
    print 'reduced chi2:', reducedChi2
    print 'pvalue:',pvalue
    print 'significance:',nSigma(1-pvalue),'sigmas'
    print ''

    fig,ax = plt.subplots()
    plotPulseProfile(ax,phaseBinEdges,phaseProfile,profileErrors,color='b',plotDoublePulse=False)
    ax.set_xlabel('phase')
    ax.set_ylabel('counts')
    ax.set_title('J0337 Pulse Profile - {}'.format(label))

if __name__=='__main__':


    paramFile = 'j0337.dict'
    params = readDict()
    params.read_from_file(paramFile)
    obsSequences = [params[key] for key in sorted(params.keys()) if key.startswith('obsSequence')]
    sunsetDates = [params[key] for key in sorted(params.keys()) if key.startswith('sunsetDate')]


#    sunsetDates.append(params['sunsetDate0'])
#    obsSequences.append(params['obsSequence0'])
#
#    sunsetDates.append(params['sunsetDate1'])
#    obsSequences.append(params['obsSequence1'])
#
#    sunsetDates.append(params['sunsetDate2'])
#    obsSequences.append(params['obsSequence2'])
    tstampsToUse = np.concatenate(obsSequences)

    nPhaseBins = 150
    wvlStart = 4000 #angstrom
    wvlEnd = 5000 #angstrom
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

    dateMasks = []
    phaseProfilesByDate = []
    for iDate,obsSeq in enumerate(obsSequences):
        dateMask = np.array([tstamp in obsSeq for tstamp in tstamps])
        phaseProfilesOnDate  = phaseProfiles[dateMask]
        combineProfiles(phaseBinEdges,phaseProfilesOnDate,label=sunsetDates[iDate])
        dateMasks.append(dateMask)
    combineProfiles(phaseBinEdges,phaseProfiles,label='total')
    
    #print '20140924',tstamps[0:17]
    #print '20140925',tstamps[17:47]
    #print '20141021',tstamps[47:]
    #profileErrors = dataDict['profileErrors']
    pop = PopUp(showMe=False)
    pop.plotArray(phaseProfiles,aspect=2.)
    pop.axes.set_yticks(np.arange(np.shape(phaseProfiles)[0]))
    pop.axes.tick_params(axis='both', which='major', labelsize=7)
    pop.axes.set_yticklabels(tstamps)
    pop.show()

    #plotArray(phaseProfiles)

#    combineProfiles(phaseBinEdges,phaseProfiles,label='total')
#    combineProfiles(phaseBinEdges,phaseProfiles[0:15],label='20140924')
#    combineProfiles(phaseBinEdges,phaseProfiles[17:47],label='20140925')
#    combineProfiles(phaseBinEdges,phaseProfiles[47:],label='20141021')
    timeProfile = np.mean(phaseProfiles,axis=1)
    fig,(ax,ax2,ax3,ax4) = plt.subplots(4,1,sharex='col')
    ax.plot(apertureRadiusList)
    ax2.plot([psfFit['flux'] for psfFit in psfFits])
    ax3.plot([psfFit['parameters'][0] for psfFit in psfFits])
    ax4.plot(timeProfile)
    for iTime,tstamp in enumerate(tstamps):
        print iTime,tstamp


    print 'done'
    plt.show()
