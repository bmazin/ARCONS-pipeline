from matplotlib import rcParams, rc
import scipy.interpolate

# common setup for matplotlib
params = {'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 13,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 13,
          #'text.usetex': True}
          'font.family':'sans-serif',
           #free font similar to Helvetica
          'font.sans-serif':'FreeSans'}

# use of Sans Serif also in math mode
rc('text.latex', preamble='\usepackage{sfmath}')

rcParams.update(params)
def cm2inch(cm):
    """Centimeters to inches"""
    return cm *0.393701

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage
import matplotlib
import mpfit

phaseShift = 1.-0.677001953125
path = '/Scratch/dataProcessing/crabData2/'
profileData = np.load(path+'fullProfile.npz')
avgProfiles=profileData['avgProfiles']
phaseBinEdges = profileData['phaseBinEdges']+phaseShift
wvlBandEdges = profileData['wvlBandEdges']
wvlBandCenters = wvlBandEdges[0:-1]+np.diff(wvlBandEdges)/2.

nOldBinsInNewBin = 66
firstAfterConvolve = nOldBinsInNewBin//2
rebinnedEdges = phaseBinEdges[::nOldBinsInNewBin]
rebinnedProfiles=scipy.ndimage.filters.convolve1d(avgProfiles,np.ones(nOldBinsInNewBin),mode='constant',axis=1)[:,firstAfterConvolve::nOldBinsInNewBin]
fig = plt.figure()
ax = fig.add_subplot(111)
print 'original shape',np.shape(avgProfiles)
print 'rebinned shape',np.shape(rebinnedProfiles)
profile = rebinnedProfiles[0]
#profile = profile / np.max(profile)
binCenters = rebinnedEdges[0:-1]+np.diff(rebinnedEdges)/2.
ax.plot(binCenters,profile,color='r')
print np.argmax(profile)

maxIdx = np.argmax(profile)
peakPhase = binCenters[maxIdx]
errorPeakPhase = binCenters[maxIdx+1]-binCenters[maxIdx]

print 'optical main pulse peak'
print peakPhase,'+/-',errorPeakPhase,' in phase'

period = 0.03362310922e6 #us

peakTime = (1.-peakPhase)*period
errorPeakTime = errorPeakPhase*period
print peakTime,'+/-',errorPeakTime,' us'

radioProfile = np.loadtxt(path+'radio/RadioProfile_LyneDM_TZRCorrect_withGUPPIdelay.txt',skiprows=1,usecols=[3])
nRadioPhaseBins = len(radioProfile)
radioProfilePhaseBins = (1.*np.arange(nRadioPhaseBins)+.5)/nRadioPhaseBins
radioProfilePhaseBins+=phaseShift
ax2 = ax.twinx()
ax2.plot(radioProfilePhaseBins,radioProfile,c=(.4,.5,.8),label='Radio Pulse')

maxRadioIdx = np.argmax(radioProfile)
radioPeakPhase = radioProfilePhaseBins[maxRadioIdx]
errorRadioPeakPhase = radioProfilePhaseBins[maxRadioIdx+1]-radioProfilePhaseBins[maxRadioIdx]

print 'radio main pulse peak'
print radioPeakPhase,'+/-',errorRadioPeakPhase, ' in phase'
print (1.-radioPeakPhase)*period,'+/-',errorRadioPeakPhase*period, ' us'

peakPhaseLag = radioPeakPhase-peakPhase

radioTimingError = 10. #us
radioTimingPhase = radioTimingError/period

errorPeakPhaseLag = np.sqrt(errorPeakPhase**2+errorRadioPeakPhase**2+radioTimingPhase**2)
peakLag = peakPhaseLag*period
errorPeakLag = errorPeakPhaseLag*period
print 'optical-radio lag'
print peakPhaseLag,'+/-',errorPeakPhaseLag, ' in phase'
print peakLag,'+/-',errorPeakLag,' us'


plt.show()



