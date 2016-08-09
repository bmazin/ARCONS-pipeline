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

def extrap1d(interpolator):
    """Taken from http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-an-extrapolated-result-beyond-the-input-range"""
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike


def inverseCdf(pdf,queryPoints):
    """Uses interpolation to invert the cdf for a discrete pdf"""
    nBins = len(pdf)
    delta = 1./(nBins-1)
    delta = queryPoints[1]-queryPoints[0]
    pdf = 1.*pdf / (delta*np.sum(pdf))
    cdf = np.cumsum(pdf)*delta
    #x = np.linspace(0.,1.,nBins)
    x = queryPoints
    #add a zero point so the valid interpolation range extends to 0
    x = np.append(0,x)
    cdf = np.append(0,cdf)
    #invCdf = scipy.interpolate.interp1d(y=x,x=cdf,bounds_error=False,fill_value=0.)
    invCdf = scipy.interpolate.interp1d(y=x,x=cdf)
    return invCdf

def inverseTransformSampler(pdf,queryPoints):
    """Returns a function that returns a given number of samples from a discrete pdf"""
    invCdf = inverseCdf(pdf,queryPoints)
    def sampler(nSamples):
        uniformProbs = np.random.random(nSamples)
        return invCdf(uniformProbs)
    return sampler

def crabSampler():
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
    return inverseTransformSampler(profile,rebinnedEdges[1:])

def sineSampler():
    #make a sinusoidal profile
    delta = 0.001
    queryPoints = np.linspace(delta,1,1/delta)
    print queryPoints[-1]
    sineProfile = np.sin(2.*np.pi*queryPoints)
    sineProfile -= np.min(sineProfile)
    return inverseTransformSampler(sineProfile,queryPoints)
    
if __name__=='__main__':
    sineSamp = sineSampler()
    
    #sineSampler = inverseTransformSampler(sineProfile,queryPoints)

    nPhotons = 1000000
    pl = sineSamp(nPhotons)
    hist,bins = np.histogram(pl,bins=500,density=True)

    crabSamp = crabSampler()
    crabPl = crabSamp(nPhotons)
    hist2,bins2 = np.histogram(crabPl,bins=500,density=True)

    fig,ax = plt.subplots(1,1)
    ax.plot(bins[0:-1],hist,color='r')
    ax.plot(bins2[0:-1],hist2,color='b')
    plt.show()
    
