#Filename:  injectPulsarKuiper.py
#Author:    Matt Strader
#Date:      May 19, 2015
#
#This script opens a list of observed photon phases, then adds fake phases
#with either a sinusoidal phase profile or a crab pulsar pulse profile.
#It adds fake phases from these profiles until the list of reak and fake phases
#shows more than N sigma signficance on a kuiper test for periodicity
import numpy as np
import tables
import numexpr
import matplotlib.pyplot as plt

from kuiper.kuiper import kuiper,kuiper_FPP
from kuiper.htest import h_test,h_fpp
from combinePulseProfiles import nSigma,plotPulseProfile
from histMetrics import kuiperFpp,hTestFpp
from fakeCrabList import fakeCrabPhases,fakeSinePhases

path = '/Scratch/dataProcessing/J0337/masterPhotons2.h5'

wvlStart = 4000.
wvlEnd = 5500.
photFile = tables.openFile(path,'r')
photTable = photFile.root.photons.photTable
i = 0
print 'cut wavelengths to range ({},{})'.format(wvlStart,wvlEnd)


phases = photTable.readWhere('(wvlStart < wavelength) & (wavelength < wvlEnd)')['phase']
photFile.close()

nPhotons = len(phases)
print nPhotons,'real photons read'

D,pval = kuiper(phases)
print 'kuiper test'
print 'D,fpp:',D,pval
print nSigma(1-pval),'sigmas'

nSigmaThreshold = 5.
nPhaseBins = 20
phaseBinEdges = np.linspace(0.,1.,nPhaseBins+1)

nSigmaSig = 0.
nFakePhotons = 0
while nSigmaSig < nSigmaThreshold:
    fakePhases = fakeSinePhases(nFakePhotons)
    groupedPhases = np.append(phases,fakePhases)

    D,pval = kuiper(groupedPhases)
    nSigmaSig = nSigma(1-pval)
    print '{} fake photons, D:{}, pval:{}, sig:{:.1f}'.format(nFakePhotons,D,pval,nSigmaSig)
    if nSigmaSig < nSigmaThreshold:
        nFakePhotons += 1000

nFakePhotonsKuiper = nFakePhotons
fakePhotonsKuiper = np.array(fakePhases)
print '{} photons to reach {} sigmas in kuiper test'.format(nFakePhotons,nSigmaThreshold)
totalProfile,_ = np.histogram(groupedPhases,bins=phaseBinEdges)
totalProfileErrors = np.sqrt(totalProfile)
profile,_ = np.histogram(phases,bins=phaseBinEdges)
profileErrors = np.sqrt(profile)

fig,ax = plt.subplots(1,1)
plotPulseProfile(ax,phaseBinEdges,profile,profileErrors=profileErrors,color='k',plotDoublePulse=False)
plotPulseProfile(ax,phaseBinEdges,totalProfile,profileErrors=totalProfileErrors,color='r',plotDoublePulse=False)
ax.set_title('Injected Sine Profile to {:.0f} sigmas - Kuiper Test'.format(nSigmaThreshold))
ax.set_xlabel('phase')
ax.set_ylabel('counts')

