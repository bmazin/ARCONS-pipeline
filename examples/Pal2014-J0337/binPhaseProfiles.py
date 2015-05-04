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

def foldPl(plPath,centroidRaStr='03:37:43.826',centroidDecStr='14:15:14.828',apertureRadius=1.,nPhaseBins=150,wvlRange=(3000,13000)):
    print 'folding',plPath
    centroidRa = float(ephem.hours(centroidRaStr))
    centroidDec = float(ephem.degrees(centroidDecStr))
    pl = PhotList(plPath)
    apertureRadius *= RADIANS_PER_ARCSEC
    wvlStart,wvlEnd = wvlRange
    print 'reading photons'
    photons = pl.photTable.readWhere('(sqrt((ra-centroidRa)**2. + (dec-centroidDec)**2.) < apertureRadius) & (wvlStart < wavelength) & (wavelength < wvlEnd)')
    pixelImageDict = makePixelImage(photons)
    
    print 'binning phases'
    del pl

    phases = photons['pulsePhase']
    phaseProfile,phaseBinEdges = np.histogram(phases,bins=nPhaseBins,range=(0.,1.))
    profileErrors = np.sqrt(phaseProfile)
    print 'done folding'
    del photons,phases
    return {'phaseProfile':phaseProfile,'phaseBinEdges':phaseBinEdges,'profileErrors':profileErrors,'pixelImage':pixelImageDict['image'],'pixelImageXRange':pixelImageDict['xRange'],'pixelImageYRange':pixelImageDict['yRange']}


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
    wvlStart = 4000 #angstrom
    wvlEnd = 5000 #angstrom
    wvlRange = (wvlStart,wvlEnd)
    centroidRaStr='03:37:43.826'
    centroidDecStr='14:15:14.828'


#    apertureData = np.loadtxt('smoothApertList.txt',dtype={'names':('obsTimestamp','oldApert','apertureRadius'),'formats':('S80','f8','f8')},delimiter='\t')
#    apertureRadiusDict = {eachItem['obsTimestamp']:eachItem['apertureRadius'] for eachItem in apertureData}
#    apertureList = [apertureRadiusDict[timestamp] for timestamp in obsTimestamps]

    foldPlBins = functools.partial(foldPl,nPhaseBins=nPhaseBins,wvlRange=wvlRange,centroidRaStr=centroidRaStr,centroidDecStr=centroidDecStr)
    def foldPlOneArg(foldArg):
        plPath,apertureRadius = foldArg
        return foldPlBins(plPath,apertureRadius=apertureRadius)

    bUseOptimalApert = True
    if bUseOptimalApert:
        dataPathApert = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_{:.1f}arcsecAperture_sigma.npz'.format(nPhaseBins,3000,8000,3.)
        apertureDataDict = np.load(dataPathApert)
        apertureList = apertureDataDict['smoothOptimalApertureRadii']
        tstamps = apertureDataDict['obsTimestamps']
        psfDicts = apertureDataDict['psfDicts']
        tstampMask = np.array([tstamp in obsTimestamps for tstamp in tstamps])
        apertureList = apertureList[tstampMask]
        tstamps = tstamps[tstampMask]
        psfDicts = psfDicts[tstampMask]
    else:
        apertureRadius=3.
        apertureList = np.ones(len(obsTstamps))*apertureRadius
    foldArgs = zip(plPaths,apertureList)

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
    outDicts = pool.map(foldPlOneArg,foldArgs)

    phaseProfiles = np.array([d['phaseProfile'] for d in outDicts])
    profileErrors = np.array([d['profileErrors'] for d in outDicts])
    pixelImages = np.array([d['pixelImage'] for d in outDicts])
    phaseBinEdges = outDicts[0]['phaseBinEdges']
    if bUseOptimalApert:
        savePath = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_optimalAperture.npz'.format(nPhaseBins,wvlStart,wvlEnd)
    else:
        savePath = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_{:.1f}arcsecAperture_sigma.npz'.format(nPhaseBins,wvlStart,wvlEnd,3.)
    np.savez(savePath,phaseProfiles=phaseProfiles,profileErrors=profileErrors,phaseBinEdges=phaseBinEdges,obsTimestamps=obsTimestamps,pixelImages=pixelImages,apertureRadiusList=apertureList,psfFits=psfDicts)
    
    fig,ax = plt.subplots()
    plotPulseProfile(ax,phaseBinEdges,phaseProfiles[0])
    print 'done'
    plt.show()
    
