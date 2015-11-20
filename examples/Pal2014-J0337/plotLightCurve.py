import numpy as np
import matplotlib.pyplot as plt
import figureHeader

def plotPulseProfile(phaseBinEdges,pulseProfile,profileErrors=None,plotDoublePulse=True,ax=None,**kwargs):
    label = kwargs.pop('label','')
    if plotDoublePulse:
        doublePhaseBinEdges = np.concatenate([phaseBinEdges,phaseBinEdges[1:]+1.])
        doubleSteppedPulseProfile = np.concatenate([pulseProfile,pulseProfile,[pulseProfile[-1]]])
        ax.plot(doublePhaseBinEdges,doubleSteppedPulseProfile,drawstyle='steps-post',label=label,**kwargs)
        if not (profileErrors is None):
            doublePulseProfile = np.concatenate([pulseProfile,pulseProfile])
            doubleProfileErrors = np.concatenate([profileErrors,profileErrors])
            doubleBinCenters = doublePhaseBinEdges[0:-1]+np.diff(doublePhaseBinEdges)/2.
            ax.errorbar(doubleBinCenters,doublePulseProfile,yerr=doubleProfileErrors,linestyle='',**kwargs)
    else:
        steppedPulseProfile = np.concatenate([pulseProfile,[pulseProfile[-1]]])
        ax.plot(phaseBinEdges,steppedPulseProfile,drawstyle='steps-post',label=label,**kwargs)
        if not (profileErrors is None):
            binCenters = phaseBinEdges[0:-1]+np.diff(phaseBinEdges)/2.
            ax.errorbar(binCenters,pulseProfile,yerr=profileErrors,linestyle='',**kwargs)


lcData = np.load('lightcurvePlot.npz')
phaseBinEdges = lcData['phaseBinEdges']
phaseProfile = lcData['phaseProfile']
profileErrors = lcData['profileErrors']
fig,ax = plt.subplots()
plotPulseProfile(phaseBinEdges,phaseProfile,profileErrors,color='k',plotDoublePulse=False,ax=ax,linewidth=1.2)
ax.set_xlabel('phase')
ax.set_ylabel('counts')

fig.savefig('lightcurve.eps')
plt.show()
