#Filename:  searchWvlCuts.py
#Author:    Matt Strader
#
#This script opens a list of photon phases and wavelengths and perferms a kuiper test on a series
# of possible wavelength cuts, to find the optimal cuts.  The metric should be compared with the results
# from running the kuiper test on the same number of phases, generated randomly.
import numpy as np
import tables
import numexpr
import functools

from kuiper.kuiper import kuiper
from kuiper.htest import h_test
from combinePulseProfiles import nSigma

path = '/Scratch/dataProcessing/J0337/masterPhotons2.h5'

def issorted(x):
    """Check if x is sorted"""
    return (np.diff(x) >= 0).all()

kuiperSorted = functools.partial(kuiper,assumeSorted=True)

wvlLimits = np.arange(3500.,11200.,200.)[0:8]
wvlStart = wvlLimits[0]
wvlEnd = wvlLimits[-1]
nWvlPairs = len(wvlLimits)**2

nTrials = 1000

metrics = np.recarray(shape=(nTrials,nWvlPairs),dtype=[('wvlStart',np.double),('wvlEnd',np.double),('D',np.double),('fpp',np.double)])
metricImg = np.zeros((nTrials,len(wvlLimits),len(wvlLimits)))
fppImg = np.zeros((nTrials,len(wvlLimits),len(wvlLimits)))
sigmaImg = np.zeros((nTrials,len(wvlLimits),len(wvlLimits)))
photFile = tables.openFile(path,'r')
photTable = photFile.root.photons.photTable
allPhotons = photTable.read()
print 'read {} photons'.format(len(allPhotons))
allWavelengths = allPhotons['wavelength']
photFile.close()

mask = numexpr.evaluate('(wvlStart < allWavelengths) & (allWavelengths < wvlEnd)')
allPhotons = allPhotons[mask]
allWavelengths = allPhotons['wavelength']
print 'initial cut to {} photons'.format(len(allPhotons))

for iTrial in np.arange(nTrials):
    print 'Trial ',iTrial
    print '***************************************************'
    if iTrial != 0:
        allPhotons['phase'] = np.random.random(len(allPhotons))

    sortIndices = np.argsort(allPhotons['phase'])
    allPhotons = allPhotons[sortIndices]
    allWavelengths = allPhotons['wavelength']
    print 'sorted phases'

    iPair = 0
    print 'sorted:',issorted(allPhotons['phase'])
    for iWvlStart,wvlStart in enumerate(wvlLimits[:-1]):
        for iWvlEnd,wvlEnd in enumerate(wvlLimits[iWvlStart+1:]):
            print 'cut wavelengths to range ({},{})'.format(wvlStart,wvlEnd)
            mask = numexpr.evaluate('(wvlStart < allWavelengths) & (allWavelengths < wvlEnd)')
            phases = allPhotons['phase'][mask]

            nPhotons = len(phases)
            print nPhotons,'photons in range'

            D,pval = kuiper(phases,assumeSorted=True)
            sig = nSigma(1-pval)
            print 'kuiper test'
            print 'D,fpp:',D,pval
            print sig,'sigmas'

            metrics[iTrial,iPair]['wvlStart'] = wvlStart
            metrics[iTrial,iPair]['wvlEnd'] = wvlEnd
            metrics[iTrial,iPair]['D'] = D
            metrics[iTrial,iPair]['fpp'] = pval
            metricImg[iTrial,iWvlStart,iWvlEnd+iWvlStart+1] = D
            fppImg[iTrial,iWvlStart,iWvlEnd+iWvlStart+1] = pval
            sigmaImg[iTrial,iWvlStart,iWvlEnd+iWvlStart+1] = sig

            del phases,mask
            iPair+=1

np.savez('/Scratch/dataProcessing/J0337/wvlCutTrials.npz',metrics=metrics,wvlLimits=wvlLimits,metricImg=metricImg,fppImg=fppImg,sigmaImg=sigmaImg)

