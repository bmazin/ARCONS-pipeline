#Filename:  fakeCrabList.py
#Author:    Matt Strader
#Date: May 19, 2015
#
#Provides a function to simulate photon phases with either a crab pulse profile
#or a sinusoidal pulse profile

import scipy.interpolate
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage
import matplotlib
import mpfit
import scipy.interpolate

phaseShift = 1.-0.677001953125
path = '/Scratch/dataProcessing/crabData2/'
profileData = np.load(path+'fullProfile1us.npz')
avgProfiles=profileData['avgProfiles']
phaseBinEdges = profileData['phaseBinEdges']#+phaseShift
wvlBandEdges = profileData['wvlBandEdges']
wvlBandCenters = wvlBandEdges[0:-1]+np.diff(wvlBandEdges)/2.


nOldBinsInNewBin = 66
firstAfterConvolve = nOldBinsInNewBin//2
rebinnedEdges = phaseBinEdges[::nOldBinsInNewBin]
rebinnedProfiles=scipy.ndimage.filters.convolve1d(avgProfiles,np.ones(nOldBinsInNewBin),mode='constant',axis=1)[:,firstAfterConvolve::nOldBinsInNewBin]
profile = rebinnedProfiles[0]
#profile = profile / np.max(profile)
profile -= np.min(profile)
profile = 1.*profile / np.sum(profile)

binCenters = rebinnedEdges[0:-1]+np.diff(rebinnedEdges)/2.

cdf = np.cumsum(profile)
nBins = len(cdf)
probs = np.linspace(0.,1.,nBins)
# y = cdf(probs)
# probs = invCdf(y)

invCdf = scipy.interpolate.interp1d(y=probs,x=cdf,bounds_error=False,fill_value=0.)

#Now for a sinusoidal profile
queryPoints = np.linspace(0.,1,10000)
sineProfile = np.sin(2.*np.pi*queryPoints)
sineProfile -= np.min(sineProfile)
sineProfile = sineProfile / np.sum(sineProfile)
sineCdf = np.cumsum(sineProfile)
nBins = len(sineCdf)
probs = np.linspace(0.,1.,nBins)
# y = cdf(probs)
# probs = invCdf(y)
sineInvCdf = scipy.interpolate.interp1d(y=probs,x=sineCdf,bounds_error=False,fill_value=0.)

def fakeCrabPhases(nPhotons=100):
    """Generates a list of phases with a profile resembling the Crab Pulsar"""
    uniformProbs = np.random.random(nPhotons)
    return invCdf(uniformProbs)

def fakeSinePhases(nPhotons=100):
    """Generates a list of phases with a sinusoidal phase profile"""
    uniformProbs = np.random.random(nPhotons)
    return sineInvCdf(uniformProbs)
    
if __name__=='__main__':
    nPhotons = 1000000
    fig,ax = plt.subplots(1,1)
    pl = fakeCrabPhases(nPhotons)
    hist,bins = np.histogram(pl,bins=500)
    hist = 1.*hist / np.sum(hist)
    ax.plot(bins[0:-1],hist)
    ax.plot(binCenters,profile)

    fig,ax = plt.subplots(1,1)
    pl = fakeSinePhases(nPhotons)
    hist,bins = np.histogram(pl,bins=500)
    hist = 1.*hist / np.sum(hist)
    ax.plot(bins[0:-1],hist)
    ax.plot(queryPoints,sineProfile)
    plt.show()
    
