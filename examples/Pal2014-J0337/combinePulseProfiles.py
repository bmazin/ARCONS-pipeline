import numpy as np
import scipy.stats,scipy.special
import matplotlib.pyplot as plt
from util.popup import plotArray
import scipy.ndimage
import os

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

    #totalProfileErrors = np.std(phaseProfiles,axis=0)
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
    nPhaseBins = 150
    wvlStart = 3000 #angstrom
    wvlEnd = 8000 #angstrom
    apertureRadius=.5#arcsec
    dataPath = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_{}arcsecAperture.npz'.format(nPhaseBins,wvlStart,wvlEnd,apertureRadius)
    print os.path.basename(dataPath)

    dataDict = np.load(dataPath)
    phaseProfiles = np.array(dataDict['phaseProfiles'],dtype=np.double)
    phaseBinEdges = dataDict['phaseBinEdges']
    tstamps = dataDict['obsTimestamps']
    #print '20140924',tstamps[0:17]
    #print '20140925',tstamps[17:47]
    #print '20141021',tstamps[47:]
    #profileErrors = dataDict['profileErrors']
    plotArray(phaseProfiles)

    combineProfiles(phaseBinEdges,phaseProfiles,label='total')
    combineProfiles(phaseBinEdges,phaseProfiles[0:15],label='20140924')
    combineProfiles(phaseBinEdges,phaseProfiles[17:47],label='20140925')
    combineProfiles(phaseBinEdges,phaseProfiles[47:],label='20141021')

    print 'done'
    plt.show()
