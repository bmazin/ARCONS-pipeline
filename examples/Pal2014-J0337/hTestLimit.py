#Filename:  hTestLimit.py
#Author:    Matt Strader
#
#This script opens a list of observed photon phases, 
import numpy as np
import tables
import numexpr
import matplotlib.pyplot as plt
import multiprocessing
import functools
import time

from kuiper.kuiper import kuiper,kuiper_FPP
from kuiper.htest import h_test,h_fpp,h_test2
from combinePulseProfiles import nSigma,plotPulseProfile
from histMetrics import kuiperFpp,hTestFpp
from inverseTransformSampling import inverseTransformSampler

def hTestTrial(iTrial,nPhotons,photonPulseFraction,pulseModel,pulseModelQueryPoints):
    np.random.seed(int((time.time()+iTrial)*1e6))
    modelSampler = inverseTransformSampler(pdf=pulseModel,queryPoints=pulseModelQueryPoints)

    nPulsePhotons = int(np.floor(photonPulseFraction*nPhotons))
    nBackgroundPhotons = int(np.ceil((1.-photonPulseFraction) * nPhotons))

    simPulsePhotons = modelSampler(nPulsePhotons)
    #background photons come from a uniform distribution
    simBackgroundPhotons = np.random.random(nBackgroundPhotons)
    simPhases = np.append(simPulsePhotons,simBackgroundPhotons)

    simHDict = h_test2(simPhases)
    simH,simM,simPval,simFourierCoeffs = simHDict['H'],simHDict['M'],simHDict['fpp'],simHDict['cs']

    print '{} - H,M,fpp,sig:'.format(iTrial),simH,simM,simPval
    return {'H':simH,'M':simM,'fpp':simPval}


path = '/Scratch/dataProcessing/J0337/masterPhotons2.h5'
trialsPath = '/Scratch/dataProcessing/J0337/randomTests4000.0-5500.0.npz'

wvlStart = 4000.
wvlEnd = 5500.
bLoadFromPl = False
hTestPath = '/Scratch/dataProcessing/J0337/hTestResults_{}-{}.npz'.format(wvlStart,wvlEnd)
nPhaseBins = 20
phaseBinEdges = np.linspace(0.,1.,nPhaseBins+1)

if bLoadFromPl:
    photFile = tables.openFile(path,'r')
    photTable = photFile.root.photons.photTable
    i = 0
    phases = photTable.readWhere('(wvlStart < wavelength) & (wavelength < wvlEnd)')['phase']
    photFile.close()
    print 'cut wavelengths to range ({},{})'.format(wvlStart,wvlEnd)

    simHs = np.load(trialsPath)['hMetrics']
    bUseSimulationFpp=False

    nPhotons = len(phases)
    print nPhotons,'real photons read'

    profile,_ = np.histogram(phases,bins=phaseBinEdges)
    profileErrors = np.sqrt(profile)

    profileErrors = profileErrors 
    profile = 1.0*profile 
    fig,ax = plt.subplots(1,1)
    plotPulseProfile(ax,phaseBinEdges,profile,profileErrors=profileErrors,color='k',plotDoublePulse=False,label='observed')
    #ax.plot(pulseModelQueryPoints,meanLevel*lightCurveModel,color='r')
    ax.set_ylabel('counts')
    ax.set_xlabel('phase')
    ax.set_title('Observed Folded Light Curve 400-550 nm')

    #Do H-test
    hDict = h_test2(phases)
    H,M,pval,fourierCoeffs = hDict['H'],hDict['M'],hDict['fpp'],hDict['cs']

    print 'h-test on real data'
    print 'H,M,fpp:',H,M,pval
    print nSigma(1-pval),'sigmas'

    #h_test2 calculates all fourierCoeffs out to 20, but for the fourier model, we only want the ones out to order M, which optimizes the Zm^2 metric
    truncatedFourierCoeffs = fourierCoeffs[0:M]
    print 'fourier coeffs:',truncatedFourierCoeffs

    #for the model, we want the negative modes as well as positve, so add them
    modelFourierCoeffs = np.concatenate([truncatedFourierCoeffs[::-1],[1.],np.conj(truncatedFourierCoeffs)])
    #make array of mode numbers
    modes = np.arange(-len(truncatedFourierCoeffs),len(truncatedFourierCoeffs)+1)

    #save so next time we can set bLoadFromPl=False
    #np.savez(hTestPath,H=H,M=M,pval=pval,fourierCoeffs=fourierCoeffs,nPhotons=nPhotons,wvlRange=(wvlStart,wvlEnd),modelFourierCoeffs=modelFourierCoeffs,modes=modes)

else:
    #Load values from previous run, when we had bLoadFromPl=True
    hTestDict = np.load(hTestPath)
    H,M,pval,fourierCoeffs,nPhotons,modelFourierCoeffs,modes = hTestDict['H'],hTestDict['M'],hTestDict['pval'],hTestDict['fourierCoeffs'],hTestDict['nPhotons'],hTestDict['modelFourierCoeffs'],hTestDict['modes']
    print 'h-test on real data'
    print 'H,M,fpp:',H,M,pval
    print nSigma(1-pval),'sigmas'
    
#make as set of x points for the pulse model we'll make
#Do NOT include x=0, or the inverted function will have a jump that causes an excess of samples
#at phase=0
nSmoothPlotPoints=1000
pulseModelQueryPoints = np.linspace(1./nSmoothPlotPoints,1,nSmoothPlotPoints)

def modelProfile(thetas):
    return np.sum( modelFourierCoeffs * np.exp(2.j*np.pi*modes*thetas[:,np.newaxis]),axis=1)
lightCurveModel = np.abs(modelProfile(pulseModelQueryPoints))
#for this test we only want the model to be the pulsed component.  We will add a DC offset later
pulseModel = lightCurveModel - np.min(lightCurveModel)

#initialPhotonPulseFraction = 1.*np.sum(pulseModel) / np.sum(lightCurveModel)
photonPulseFraction=15400./nPhotons #skip to previously determined answer
print 'photon fraction',photonPulseFraction

#get samples with distribution of the modelProfile
#modelSampler = inverseTransformSampler(pdf=lightCurveModel,queryPoints=pulseModelQueryPoints)
modelSampler = inverseTransformSampler(pdf=pulseModel,queryPoints=pulseModelQueryPoints)

nTrials = 1
#for each trial run the h test on a set of photon phases with our model profile, and with the pulse fraction specified
#we want to make a distribution of H values for this pulse fraction, model, and number of photons

#make a function that only takes the trial number (as an identifier)
mappableHTestTrial = functools.partial(hTestTrial,pulseModel=pulseModel,
        pulseModelQueryPoints=pulseModelQueryPoints,nPhotons=nPhotons,
        photonPulseFraction=photonPulseFraction)
pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-3)#leave a few processors for other people
outDicts = pool.map(mappableHTestTrial,np.arange(nTrials))

simHs = np.array([out['H'] for out in outDicts])
simPvals = np.array([out['fpp'] for out in outDicts])
#save the resulting list of H vals
np.savez('sim-h-{}.npz'.format(nTrials),simHs=simHs,simPvals=simPvals,pval=pval,H=H,photonPulseFraction=photonPulseFraction,nPhotons=nPhotons)

#make a model profile once more for a plot
modelSampler = inverseTransformSampler(pdf=pulseModel,queryPoints=pulseModelQueryPoints)

nPulsePhotons = int(np.floor(photonPulseFraction*nPhotons))
nBackgroundPhotons = int(np.ceil((1.-photonPulseFraction) * nPhotons))

simPulsePhotons = modelSampler(nPulsePhotons)
#background photons come from a uniform distribution
simBackgroundPhotons = np.random.random(nBackgroundPhotons)
#put them together for the full profile
simPhases = np.append(simPulsePhotons,simBackgroundPhotons)

#run the H test on the generated phases
simHDict = h_test2(simPhases)
simH,simM,simPval,simFourierCoeffs = simHDict['H'],simHDict['M'],simHDict['fpp'],simHDict['cs']

print 'h-test on simulated samples'
print 'H,M,fpp:',simH,simM,simPval
print nSigma(1-simPval),'sigmas'


simProfile,_ = np.histogram(simPhases,bins=phaseBinEdges)
simProfileErrors = np.sqrt(simProfile)
meanLevel = np.mean(simProfile)
print 'mean level',meanLevel

#plotPulseProfile(ax,phaseBinEdges,profile,profileErrors=profileErrors,color='b',plotDoublePulse=False,label='sim')
#
#plt.show()

print 'fraction of photons in pulse',photonPulseFraction
print '{} trials'.format(len(simHs))
print 'observed H',H

frac = 1.*np.sum(simPvals<pval)/len(simPvals)
print 'fraction below observed fpp:',frac
#hHist,hBinEdges = np.histogram(simHs,bins=100,density=True)
fppHist,fppBinEdges = np.histogram(simPvals,bins=100,density=True)

fig,ax = plt.subplots(1,1)
ax.plot(fppBinEdges[0:-1],fppHist,drawstyle='steps-post',color='k')
ax.axvline(pval,color='r')
ax.set_xlabel('fpp')
ax.set_ylabel('frequency')

ax.set_title('H-test')

magG = 17.93
sineMagDiff = -2.5*np.log10(photonPulseFraction)
print 'SDSS magnitude g: {:.2f}'.format(magG)
print 'magnitude difference: {:.2f}'.format(sineMagDiff)
print 'limiting g mag: {:.2f}'.format(magG+sineMagDiff)

plt.show()

