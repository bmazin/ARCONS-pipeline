import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,FixedLocator
from util.popup import plotArray,PopUp
from combinePulseProfiles import nSigma
from kuiper.htest import h_fpp
#metrics=metrics,wvlLimits=wvlLimits,metricImg=metricImg,fppImg=fppImg,sigmaImg=sigmaImg,wvlPairs=wvlPairs,energyBinWidth=energyBinWidth)
dirPath = '/Scratch/dataProcessing/J0337/'
#globPath = os.path.join(dirPath,'chtest/chtest*')
globPath = os.path.join(dirPath,'chtestAll.txt')

metricImgs = []
fppImgs = []
sigmaImgs = []
linesPerTrial = 276
bestRealH = 12.2699
for iPath,path in enumerate(glob.glob(globPath)):

    data = np.loadtxt(path)
    if (len(data) > 0) and ((len(data) % linesPerTrial) == 0):
        wvlStartList = data[:,0]
        wvlStopList = data[:,1]
        hMetrics = data[:,2]

        wvlPairs = zip(wvlStartList,wvlStopList)
        wvlLimits = np.unique(np.append(wvlStartList,wvlStopList))
        nWvlPairs = linesPerTrial #len(wvlPairs)
        nWvlStarts = len(wvlLimits) - 1
        nWvlLimits = len(wvlLimits)

        iWvlPair = 0
        for iWvlSet in xrange(len(data)//nWvlPairs):
            metricImg = np.zeros((nWvlLimits,nWvlLimits))
            fppImg = np.zeros((nWvlLimits,nWvlLimits))
            sigmaImg = np.zeros((nWvlLimits,nWvlLimits))
            for iWvlStart,wvlStart in enumerate(wvlLimits[:-1]):
                for iWvlEnd,wvlEnd in enumerate(wvlLimits[iWvlStart+1:]):
                    assert wvlStart == wvlStartList[iWvlPair]
                    assert wvlEnd == wvlStopList[iWvlPair]

                    h = hMetrics[iWvlPair]
                    pval = h_fpp(h)
                    sig = nSigma(1-pval)
                    
                    metricImg[iWvlStart,iWvlEnd+iWvlStart+1] = h
                    fppImg[iWvlStart,iWvlEnd+iWvlStart+1] = pval
                    sigmaImg[iWvlStart,iWvlEnd+iWvlStart+1] = sig
                    iWvlPair += 1
            metricImgs.append(metricImg)
            fppImgs.append(fppImg)
            sigmaImgs.append(sigmaImg)

metricImgs = np.array(metricImgs)
fppImgs = np.array(fppImgs)
sigmaImgs = np.array(sigmaImgs)
metricImgs = metricImgs[0:1000] #trim to just 1000 trials

nTrials = len(metricImgs)

bestSig = np.max(sigmaImgs)
print '{} trials'.format(nTrials)

print 'highest:',bestSig,'sigmas'
metricsAboveBest = (metricImgs >= bestRealH)
#metricsAboveBest = (sigmaImgs >= 3.2)
trialsAboveBest = np.any(metricsAboveBest,axis=(1,2))
nTrialsAboveBest = np.sum(trialsAboveBest)

fractionAboveBest = 1.*nTrialsAboveBest/nTrials
print fractionAboveBest, 'of metrics above the real best h',bestRealH
print 'so significance of best h is'
print nSigma(1-fractionAboveBest),'sigmas'

np.savez('/Scratch/dataProcessing/J0337/cWvlTrials_{}.npz'.format(nTrials),metricImgs=metricImgs,fppImgs=fppImgs,sigmaImgs=sigmaImgs,wvlLimits=wvlLimits,nTrials=nTrials,fractionAboveBest=fractionAboveBest,bestRealh=bestRealH)
print 'done!'





