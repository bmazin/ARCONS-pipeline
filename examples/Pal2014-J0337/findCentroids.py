import numpy as np
import tables
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
import os
import hotpix.hotPixels as hp
from util.FileName import FileName
from util.popup import plotArray
from util.readDict import readDict
import astrometry.CentroidCalc as cc
import multiprocessing

def centroidObs(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess):
    obs = ObsFile(obsPath)
    print obsPath,obs.getFromHeader('exptime'),obs
    obs.loadAllCals()
#    obs.loadBestWvlCalFile()
#    obs.loadFlatCalFile(flatPath)
    obs.setWvlCutoffs(3000,11000)
#    if not os.path.exists(hotPath):
#        hp.findHotPixels(obsFile=obs,outputFileName=hotPath,display=True,fwhm=2.,boxSize=5, nSigmaHot=4.0,)
#    obs.loadHotPixCalFile(hotPath,switchOnMask=True)
    cc.centroidCalc(obs,centroidRa,centroidDec,guessTime=300,integrationTime=30,secondMaxCountsForDisplay=2000,HA_offset=haOffset,xyapprox=[xGuess,yGuess],outputFileName=centroidPath,usePsfFit=True,radiusOfSearch=8)
    print 'done centroid',centroidPath
    del obs

if __name__=='__main__':
    paramFile = 'j0337.dict'
    params = readDict()
    params.read_from_file(paramFile)
    run = params['run']
    sunsetDates = []
    flatDates = []
    obsSequences = []
    parFiles = []
    guessRows = []
    guessCols = []
    centerPixelGuesses = []

    sunsetDates.append(params['sunsetDate0'])
    flatDates.append(params['flatDate0'])
    parFiles.append(params['parFile0'])
    obsSequences.append(params['obsSequence0'])
    centerPixelGuesses.append(params['centerPixel0'])

    sunsetDates.append(params['sunsetDate1'])
    flatDates.append(params['flatDate1'])
    parFiles.append(params['parFile1'])
    obsSequences.append(params['obsSequence1'])
    centerPixelGuesses.append(params['centerPixel1'])

    sunsetDates.append(params['sunsetDate2'])
    flatDates.append(params['flatDate2'])
    parFiles.append(params['parFile2'])
    obsSequences.append(params['obsSequence2'])
    centerPixelGuesses.append(params['centerPixel2'])
    


    obsFilenames = []
    timeMaskFilenames = []
    flatFilenames = []
    centroidFilenames = []
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        sunsetDate = sunsetDates[iSeq]
        obsFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).obs() for ts in obsSequence])
        timeMaskFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).timeMask() for ts in obsSequence])
        flatFilenames.append(FileName(run=run,date=flatDates[iSeq],tstamp='').flatSoln())
        centroidFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).centroidList() for ts in obsSequence])


    centroidRa = '03:37:43.826'
    centroidDec = '14:15:14.828'
    haOffset = 150.
    print flatFilenames

    nWorkers = 7
    activeWorkers = []
    print centerPixelGuesses
    for iSeq,obsSeq in enumerate(obsFilenames):
        for iObs,obsPath in enumerate(obsSeq):
            print 'starting',obsPath
            centroidPath = centroidFilenames[iSeq][iObs]
            yGuess,xGuess = centerPixelGuesses[iSeq]
            flatPath = flatFilenames[iSeq]
            if nWorkers > 1:
                worker = multiprocessing.Process(target=centroidObs,args=(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess))
                worker.start()
                activeWorkers.append(worker)
            else:
                print 'entering centroid'
                centroidObs(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess)
                print 'exiting centroid'
            if len(activeWorkers) >= nWorkers:
                for worker in activeWorkers:
                    worker.join()
                activeWorkers = []

    
