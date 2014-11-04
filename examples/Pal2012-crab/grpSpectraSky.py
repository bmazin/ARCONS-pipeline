from matplotlib import rcParams, rc

# common setup for matplotlib
params = {'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 14,
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
from util.ObsFile import ObsFile
import tables
from astropy import units as u
from astropy import constants as c
import scipy.stats
import scipy.special

def nSigma(pvalue):
    return scipy.special.erfinv(pvalue)*np.sqrt(2.)

np.seterr(divide='ignore')

path = '/Scratch/dataProcessing/crabData2/'
nIdxToCheck = 81
nSigmaRadioCutoff = 3
nBins=250
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_skyW.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_KS.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
#dataFilePath = path+'indPulseProfiles_{}sigma_{}_{}phaseBins_KSswapStrong.h5'.format(nSigmaRadioCutoff,nIdxToCheck,nBins)
dataFilePath = path+'indPulseProfiles_{}sigma_P2_KS.h5'.format(nSigmaRadioCutoff)
dataFile = tables.openFile(dataFilePath,mode='r')
idxOffsets = dataFile.root.idxOffsets.read()
idx0 = np.searchsorted(idxOffsets,0)

counts = dataFile.root.counts.read()
wvlBinEdges = dataFile.root.wvlBinEdges.read()
fullSpectra = dataFile.root.fullSpectra.read()
skyFullSpectra = dataFile.root.skyFullSpectra.read()
peakSpectra = dataFile.root.peakSpectra.read()
skyPeakSpectra = dataFile.root.skyPeakSpectra.read()
radioIndices = dataFile.root.radioIndices.read()
radioPhaseMask = np.logical_and(radioIndices >= 1369,radioIndices <= 1394)
earlyPhaseMask = np.logical_and(radioIndices >= 1373,radioIndices <= 1378)

#dataFile.root.grpWavelengths.read()[iGRP,iPhoton]

wvlBinWidths = np.diff(wvlBinEdges)
wvlBinCenters = wvlBinEdges[0:-1]+wvlBinWidths/2.

nIdxOffsets = len(idxOffsets)
cmap = matplotlib.cm.jet
nPulsesPerIdx = np.sum(counts>0,axis=0)


dimMask = np.ones(len(counts))
idx0 = np.searchsorted(idxOffsets,0)
dimMask[counts[:,idx0]==0]=0
dimMask = (dimMask==1)

radioStrength = dataFile.root.radioMax.read()
radioStrength = radioStrength[dimMask]
radioStrengthCutoff = .175#0.175
strongMask = np.array(radioStrength >= radioStrengthCutoff)

radioPhaseMask = radioPhaseMask[dimMask][strongMask]
earlyPhaseMask = earlyPhaseMask[dimMask][strongMask]

grpWavelengths = dataFile.root.grpWavelengths.read()
grpWavelengths = np.array(grpWavelengths)
nongrpWavelengths = dataFile.root.nongrpWavelengths.read()
nongrpWavelengths = np.array(nongrpWavelengths)[strongMask]

print np.sum(radioPhaseMask)
grpWavelengthList = np.concatenate((grpWavelengths[radioPhaseMask]).tolist())
earlygrpWavelengthList = np.concatenate(grpWavelengths[earlyPhaseMask].tolist())

nongrpWavelengths = np.concatenate(nongrpWavelengths[radioPhaseMask].tolist())

ks = scipy.stats.ks_2samp(grpWavelengthList,nongrpWavelengths)
print 'KS test GRP,nonGRP',ks,nSigma(1-ks[1]),'sigmas'
ks = scipy.stats.ks_2samp(grpWavelengthList,earlygrpWavelengthList)
print 'KS test GRP,early GRP',ks,nSigma(1-ks[1]),'sigmas'
ks = scipy.stats.ks_2samp(nongrpWavelengths,earlygrpWavelengthList)
print 'KS test nonGRP,early GRP',ks,nSigma(1-ks[1]),'sigmas'
fig = plt.figure()
ax = fig.add_subplot(111)
spectraToAvg = []
peakSpectraToAvg = []
skyPeakSpectraToAvg = []

for iIdxOffset,idxOffset in enumerate(idxOffsets):
    colorFraction = (iIdxOffset+1.)/nIdxOffsets
    color=cmap(colorFraction)
    if idxOffset == 0:
        color = (0.,0.,0.)
    fullSpectrum = fullSpectra[iIdxOffset]/nPulsesPerIdx[iIdxOffset]
    peakSpectrum = peakSpectra[iIdxOffset]/nPulsesPerIdx[iIdxOffset]
    skyPeakSpectrum = skyPeakSpectra[iIdxOffset]/nPulsesPerIdx[iIdxOffset]
    if idxOffset == 0:
        fullSpectrum0 = fullSpectrum
        peakSpectrum0 = peakSpectrum
    else:
        spectraToAvg.append(fullSpectrum)
        peakSpectraToAvg.append(peakSpectrum)
    skyPeakSpectraToAvg.append(skyPeakSpectrum)
    ax.plot(wvlBinCenters,peakSpectrum/(1.*np.sum(peakSpectrum)),c=color,label='%d'%idxOffset)
ax.plot(wvlBinCenters,peakSpectrum0/(1.*np.sum(peakSpectrum0)),c='k',label='%d'%idxOffset)

fig = plt.figure()
ax = fig.add_subplot(111)
spectraToAvg = np.array(spectraToAvg)
avgSpectrum = np.mean(spectraToAvg,axis=0)
stdSpectrum = np.std(spectraToAvg,axis=0)
ax.plot(wvlBinCenters,avgSpectrum,label='surrounding pulses',color='r')
ax.errorbar(wvlBinCenters,fullSpectrum0,yerr=stdSpectrum,label='GRP-accompanied',color='k')

ax.set_xlabel('Wavelength ($\AA$)')
ax.set_ylabel('Counts per $\AA$ per pixel per pulse')

# Put a legend below current axis
ax.legend(loc='best')

peakSpectraToAvg = np.array(peakSpectraToAvg)
avgSpectrum = np.mean(peakSpectraToAvg,axis=0)
stdSpectrum = np.std(peakSpectraToAvg,axis=0)

fig = plt.figure()
ax = fig.add_subplot(111)
#avgSpectrum /= np.sum(avgSpectrum)
#stdSpectrum /= np.sum(peakSpectrum)
#peakSpectrum0 /= np.sum(peakSpectrum)
ax.plot(wvlBinCenters,avgSpectrum,label='surrounding pulses',color='r')
ax.errorbar(wvlBinCenters,peakSpectrum0,yerr=stdSpectrum,label='GRP-accompanied',color='k')

ax.set_xlabel('Wavelength ($\AA$)')
ax.set_ylabel('Counts per $\AA$ per pixel per pulse')
ax.set_title('Spectrum at Main Peak (Not Sky Subtracted)')

skyPeakSpectraToAvg = np.array(skyPeakSpectraToAvg)
avgSkySpectrum = np.mean(skyPeakSpectraToAvg,axis=0)
stdSkySpectrum = np.std(skyPeakSpectraToAvg,axis=0)
ax.plot(wvlBinCenters,avgSkySpectrum,label='sky',color='b')
ax.legend(loc='best')

fig = plt.figure()
ax = fig.add_subplot(111)
skySubAvgPeakSpectrum = avgSpectrum-avgSkySpectrum
skySubPeakSpectrum0 = peakSpectrum0-avgSkySpectrum
#Still need to adjust errors for subtraction
skySubAvgPeakSpectrum /= np.sum(skySubAvgPeakSpectrum)
stdSpectrum /= np.sum(skySubPeakSpectrum0)
skySubPeakSpectrum0 /= np.sum(skySubPeakSpectrum0)
#ax.plot(wvlBinCenters,skySubAvgPeakSpectrum,label='surrounding non-GRP-coincident',color='r',linestyle='steps-mid',linewidth=1.5)
ax.plot(wvlBinEdges,np.append(skySubAvgPeakSpectrum,skySubAvgPeakSpectrum[-1]),'r',drawstyle='steps-post',label='surrounding non-GRP-accompanied',linewidth=1.5)
ax.errorbar(wvlBinCenters,skySubPeakSpectrum0,yerr=stdSpectrum,color='k',linestyle='',linewidth=1.5)
ax.plot(wvlBinEdges,np.append(skySubPeakSpectrum0,skySubPeakSpectrum0[-1]),'k',drawstyle='steps-post',label='GRP-accompanied',linewidth=1.5)

ax.set_xlabel('Wavelength ($\AA$)')
ax.set_ylabel('Normalized Counts per Wavelength Bin')
ax.legend(loc='best')
#ax.set_title('Spectrum at Main Peak')


plt.show()


