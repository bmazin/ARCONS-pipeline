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
trialsPath = '/Scratch/dataProcessing/J0337/randomTests4000.0-5500.0.npz'

wvlStart = 4000.
wvlEnd = 5500.
photFile = tables.openFile(path,'r')
photTable = photFile.root.photons.photTable
i = 0
print 'cut wavelengths to range ({},{})'.format(wvlStart,wvlEnd)
phases = photTable.readWhere('(wvlStart < wavelength) & (wavelength < wvlEnd)')['phase']
photFile.close()

simDs = np.load(trialsPath)['dMetrics']
bUseSimulationFpp=False

nPhotons = len(phases)
print nPhotons,'real photons read'

D,pval = kuiper(phases)
simFpp = 1.*np.sum(simDs>=D) / len(simDs)

print 'kuiper test'
print 'D,fpp:',D,pval
print 'simulated fpp:',simFpp
print nSigma(1-pval),'sigmas'
print nSigma(1-simFpp),'sigmas (sim)'

nSigmaThreshold = 5.
nPhaseBins = 20
phaseBinEdges = np.linspace(0.,1.,nPhaseBins+1)
profile,_ = np.histogram(phases,bins=phaseBinEdges)
profileErrors = np.sqrt(profile)

magG = 17.93
photonIncrement = 1000

#Inject Sine first
nSigmaSig = 0.
nFakePhotons = 0
while nSigmaSig < nSigmaThreshold:
    fakePhases = fakeSinePhases(nFakePhotons)
    groupedPhases = np.append(phases,fakePhases)

    D,pval = kuiper(groupedPhases)
    if bUseSimulationFpp:
        pval = 1.*np.sum(simDs>=D) / len(simDs)
    nSigmaSig = nSigma(1-pval)
    #print '{} fake photons, D:{}, pval:{}, sig:{:.1f}'.format(nFakePhotons,D,pval,nSigmaSig)
    if nSigmaSig < nSigmaThreshold:
        nFakePhotons += photonIncrement

nFakeSinePhotons = nFakePhotons
fakeSinePhotons= np.array(fakePhases)
print '{} sine-profile photons to reach {} sigmas in kuiper test'.format(nFakePhotons,nSigmaThreshold)
sineMagDiff = -2.5*np.log10(1.*nFakePhotons/nPhotons)
print 'magnitude difference: {:.2f}'.format(sineMagDiff)
print 'limiting g mag: {:.2f}'.format(magG+sineMagDiff)

totalSineProfile,_ = np.histogram(groupedPhases,bins=phaseBinEdges)
totalSineProfileErrors = np.sqrt(totalSineProfile)

nSigmaSig = 0.
nFakePhotons = 0
while nSigmaSig < nSigmaThreshold:
    fakePhases = fakeCrabPhases(nFakePhotons)
    groupedPhases = np.append(phases,fakePhases)

    D,pval = kuiper(groupedPhases)
    if bUseSimulationFpp:
        pval = 1.*np.sum(simDs>=D) / len(simDs)
    nSigmaSig = nSigma(1-pval)
    #print '{} fake photons, D:{}, pval:{}, sig:{:.1f}'.format(nFakePhotons,D,pval,nSigmaSig)
    if nSigmaSig < nSigmaThreshold:
        nFakePhotons += photonIncrement

nFakeCrabPhotons= nFakePhotons
fakeCrabPhotons= np.array(fakePhases)
print '{} crab-profile photons to reach {} sigmas in kuiper test'.format(nFakePhotons,nSigmaThreshold)
crabMagDiff = -2.5*np.log10(1.*nFakePhotons/nPhotons)
print 'magnitude difference: {:.2f}'.format(crabMagDiff)
print 'limiting g mag: {:.2f}'.format(magG+crabMagDiff)

totalCrabProfile,_ = np.histogram(groupedPhases,bins=phaseBinEdges)
totalCrabProfileErrors = np.sqrt(totalCrabProfile)

fig,axs = plt.subplots(2,1)
plotPulseProfile(axs[0],phaseBinEdges,profile,profileErrors=profileErrors,color='k',plotDoublePulse=False,label='observed')
plotPulseProfile(axs[0],phaseBinEdges,totalSineProfile,profileErrors=totalSineProfileErrors,color='r',plotDoublePulse=False,label='injected with sine profile')

plotPulseProfile(axs[1],phaseBinEdges,profile,profileErrors=profileErrors,color='k',plotDoublePulse=False,label='observed')
plotPulseProfile(axs[1],phaseBinEdges,totalCrabProfile,profileErrors=totalCrabProfileErrors,color='r',plotDoublePulse=False,label='injected with crab profile')
axs[0].set_title('Injected Photons to {:.0f} sigmas - Kuiper Test'.format(nSigmaThreshold))
axs[0].set_ylabel('counts')
axs[0].legend(loc='best')
axs[1].set_ylabel('counts')
axs[1].set_xlabel('phase')
axs[1].legend(loc='best')

plt.show()

