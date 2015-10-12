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

wvlLimits = np.arange(3500.,11200.,200.)[0:7]
wvlStart = wvlLimits[0]
wvlEnd = wvlLimits[-1]
nWvlPairs = len(wvlLimits)**2
metrics = np.recarray(shape=(nWvlPairs,),dtype=[('wvlStart',np.double),('wvlEnd',np.double),('D',np.double),('fpp',np.double)])
metricImg = np.zeros((len(wvlLimits),len(wvlLimits)))
fppImg = np.zeros((len(wvlLimits),len(wvlLimits)))
sigmaImg = np.zeros((len(wvlLimits),len(wvlLimits)))
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

if bRandom:
    allPhotons['phase'] = np.random.random(len(allPhotons))

if ~issorted(allPhotons['phase']):
    sortIndices = np.argsort(allPhotons['phase'])
    allPhotons = allPhotons[sortIndices]
    print 'sorted phases'


i = 0
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

        metrics[i]['wvlStart'] = wvlStart
        metrics[i]['wvlEnd'] = wvlEnd
        metrics[i]['D'] = D
        metrics[i]['fpp'] = pval
        metricImg[iWvlStart,iWvlEnd+iWvlStart+1] = D
        fppImg[iWvlStart,iWvlEnd+iWvlStart+1] = pval
        sigmaImg[iWvlStart,iWvlEnd+iWvlStart+1] = sig

        del phases
        i+=1

np.savez('/Scratch/dataProcessing/J0337/wvlCutMetricList.npz',metrics=metrics,wvlLimits=wvlLimits,metricImg=metricImg,fppImg=fppImg,sigmaImg=sigmaImg)

