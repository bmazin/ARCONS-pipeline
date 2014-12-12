import numpy as np
import scipy.stats,scipy.special
import matplotlib.pyplot as plt
from util.popup import plotArray

def plotPulseProfile(ax,phaseBinEdges,pulseProfile,profileErrors=None,**kwargs):
    doublePhaseBinEdges = np.concatenate([phaseBinEdges,phaseBinEdges[1:]+1.])
    doubleSteppedPulseProfile = np.concatenate([pulseProfile,pulseProfile,[pulseProfile[-1]]])
    ax.plot(doublePhaseBinEdges,doubleSteppedPulseProfile,drawstyle='steps-post',**kwargs)
    if not (profileErrors is None):
        doublePulseProfile = np.concatenate([pulseProfile,pulseProfile])
        doubleProfileErrors = np.concatenate([profileErrors,profileErrors])
        doubleBinCenters = doublePhaseBinEdges[0:-1]+np.diff(doublePhaseBinEdges)/2.
        ax.errorbar(doubleBinCenters,doublePulseProfile,yerr=doubleProfileErrors,linestyle='',**kwargs)

def nSigma(pvalue):
    return scipy.special.erfinv(pvalue)*np.sqrt(2.)

if __name__=='__main__':
    dataPath = '/Scratch/dataProcessing/J0337/profiles201409.npz'
    dataDict = np.load(dataPath)
    phaseProfiles = np.array(dataDict['phaseProfiles'],dtype=np.double)
    phaseBinEdges = dataDict['phaseBinEdges']
    profileErrors = dataDict['profileErrors']

    print np.shape(phaseProfiles)
    #phaseProfiles = phaseProfiles[0:39]
    #phaseProfiles = phaseProfiles[0:17]
    #phaseProfiles = phaseProfiles[14:16]
    #phaseProfiles = phaseProfiles[17:39]
    plotArray(phaseProfiles)

    phaseProfile = np.sum(phaseProfiles,axis=0)

    #totalProfileErrors = np.std(phaseProfiles,axis=0)
    profileErrors = np.sqrt(phaseProfile)

    profileAvg = np.average(phaseProfile,weights=1./profileErrors**2)
    print profileAvg
    errorProfileAvg = 1./np.sqrt(np.sum(1./profileErrors**2))
    print 'weighted average(%):',profileAvg,'+/-',errorProfileAvg
    chi2=np.sum((phaseProfile-profileAvg)**2/profileErrors**2)
    dof = len(phaseProfile) - 1 #1 free parameter - the average
    reducedChi2 = 1.*chi2/dof
    pvalue=1-scipy.stats.chi2.cdf(chi2,dof)
    print 'compare to flat line'
    print 'chi2:', chi2
    print 'reduced chi2:', reducedChi2
    print 'pvalue:',pvalue
    print 'significance:',nSigma(1-pvalue),'sigmas'

    fig,ax = plt.subplots()
    plotPulseProfile(ax,phaseBinEdges,phaseProfile,profileErrors,color='b')
    ax.set_xlabel('phase')
    ax.set_ylabel('counts')
    ax.set_title('J0337 Pulse Profile')
    print 'done'
    plt.show()
