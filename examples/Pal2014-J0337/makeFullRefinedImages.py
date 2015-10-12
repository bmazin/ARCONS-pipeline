import numpy as np
import tables
import matplotlib.pyplot as plt
from util.FileName import FileName
from util.popup import plotArray,PopUp
import multiprocessing
from photonlist.photlist import writePhotonList,PhotList,xyPack
from headers.ArconsHeaders import PulsarPhotonList
from timing import photonTiming
import ephem
import functools
from util.readDict import readDict
import os

RADIANS_PER_ARCSEC = np.pi/648000.

def plotPulseProfile(ax,phaseBinEdges,pulseProfile,**kwargs):
    doublePhaseBinEdges = np.append(phaseBinEdges,phaseBinEdges[1:]+1.)
    doublePulseProfile = np.append(pulseProfile,pulseProfile)
    ax.plot(doublePhaseBinEdges,np.append(doublePulseProfile,doublePulseProfile[-1]),drawstyle='steps-post',**kwargs)

def makePixelImage(photons):
    xPix = photons['xPix']
    yPix = photons['yPix']
    xRange = np.arange(np.min(xPix),np.max(xPix)+1)
    yRange = np.arange(np.min(yPix),np.max(yPix)+1)

    xyCounts = np.bincount(photons['xyPix'],minlength=xyPack(col=xRange[-1],row=yRange[-1])+1)
    image = np.zeros((46,44))
    for iX,x in enumerate(xRange):
        for iY,y in enumerate(yRange):
            counts = xyCounts[xyPack(row=y,col=x)]
            image[y,x] = counts
    return {'image':image,'xRange':xRange,'yRange':yRange}

def extractImage(plPath,centroidRaStr='03:37:43.826',centroidDecStr='14:15:14.828',apertureRadius=1.,wvlRange=(3000,13000)):
    print 'extracting',plPath
    centroidRa = float(ephem.hours(centroidRaStr))
    centroidDec = float(ephem.degrees(centroidDecStr))
    pl = PhotList(plPath)
    apertureRadius *= RADIANS_PER_ARCSEC
    wvlStart,wvlEnd = wvlRange
    print 'reading photons'
    photons = pl.photTable.readWhere('(sqrt((ra-centroidRa)**2. + (dec-centroidDec)**2.) < apertureRadius) & (wvlStart < wavelength) & (wavelength < wvlEnd)')
    pixelImageDict = makePixelImage(photons)
    
    del pl

    del photons
    return {'pixelImage':pixelImageDict['image'],'pixelImageXRange':pixelImageDict['xRange'],'pixelImageYRange':pixelImageDict['yRange']}


if __name__=='__main__':
    savePath = '/Scratch/dataProcessing/J0337/'
    paramFile = 'j0337.dict'
    params = readDict()
    params.read_from_file(paramFile)
    run = params['run']
    sunsetDates = []
    flatDates = []
    obsSequences = []
    parFiles = []

    sunsetDates.append(params['sunsetDate0'])
    flatDates.append(params['flatDate0'])
    parFiles.append(params['parFile0'])
    obsSequences.append(params['obsSequence0'])
    
    sunsetDates.append(params['sunsetDate1'])
    flatDates.append(params['flatDate1'])
    parFiles.append(params['parFile1'])
    obsSequences.append(params['obsSequence1'])

    sunsetDates.append(params['sunsetDate2'])
    flatDates.append(params['flatDate2'])
    parFiles.append(params['parFile2'])
    obsSequences.append(params['obsSequence2'])


    plPaths = []
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        sunsetDate = sunsetDates[iSeq]
        plPaths.append([FileName(run=run,date=sunsetDate,tstamp=ts).photonList() for ts in obsSequence])

    obsTimestamps = np.concatenate(obsSequences)
    plPaths = [path for pathList in plPaths for path in pathList]

    nPhaseBins = 150 #bin it finely, we can rebin later
    wvlStart = 3000 #angstrom
    wvlEnd = 8000 #angstrom
    wvlRange = (wvlStart,wvlEnd)
    apertureRadius=6.#arcsec
    extractImages = functools.partial(extractImage,wvlRange=wvlRange,apertureRadius=apertureRadius)

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
    outDicts = pool.map(extractImages,plPaths)

    pixelImages = np.array([d['pixelImage'] for d in outDicts])

    savePath = '/Scratch/dataProcessing/J0337/plImages_2014_{}bins_{}-{}angstroms_{}arcsecAperture.npz'.format(nPhaseBins,wvlStart,wvlEnd,apertureRadius)
    np.savez(savePath,obsTimestamps=obsTimestamps,pixelImages=pixelImages)
    
    print 'done'
    plt.show()
    
