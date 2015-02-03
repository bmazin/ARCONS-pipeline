import numpy as np
import tables
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
import os
import hotpix.hotPixels as hp
from util.FileName import FileName
from util.popup import plotArray,PopUp
from util.readDict import readDict
import astrometry.CentroidCalc as cc
import multiprocessing

def centroidObs(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess,savePath,tstamp):
    obs = ObsFile(obsPath)
#    if not os.path.exists(hotPath):
#        hp.findHotPixels(obsFile=obs,outputFileName=hotPath)
    obs.loadAllCals()
#    obs.loadHotPixCalFile(hotPath,switchOnMask=True)
#    obs.loadBestWvlCalFile()
#    obs.loadFlatCalFile(flatPath)
    obs.setWvlCutoffs(3000,11000)
    obs.loadCentroidListFile(centroidPath)
    
    ctrdFile = obs.centroidListFile
    sliceTimes = ctrdFile.root.centroidlist.times.read()
    xPositions = ctrdFile.root.centroidlist.xPositions.read()
    yPositions = ctrdFile.root.centroidlist.yPositions.read()
    intTime = sliceTimes[1]-sliceTimes[0]
    
    for iTime,time in enumerate(sliceTimes):
        x = xPositions[iTime]
        y = yPositions[iTime]
        title='centroid_{}_{}s'.format(tstamp,time)
        imgDict = obs.getPixelCountImage(firstSec=time,integrationTime=intTime,weighted=True)
        imgPath=os.path.join(savePath,title+'.png')
        pop = PopUp(showMe=False)
        pop.plotArray(imgDict['image'],title=title)
        pop.axes.plot(x,y,color='g',marker='d')
        pop.fig.savefig(imgPath)
        print 'saved to',imgPath
        
    del obs

if __name__=='__main__':
    savePath = '/Scratch/dataProcessing/J0337/'
    run = 'PAL2014'
    sunsetDates = []
    flatDates = []
    obsSequences = []
    guessRows = []
    guessCols = []
    
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
    parFiles.append(params['parFile0'])
    obsSequences.append(params['obsSequence0'])
    centerPixelGuesses.append(params['centerPixel0'])

    sunsetDates.append(params['sunsetDate1'])
    parFiles.append(params['parFile1'])
    obsSequences.append(params['obsSequence1'])
    centerPixelGuesses.append(params['centerPixel1'])

    sunsetDates.append(params['sunsetDate2'])
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
        #plFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).crabList() for ts in obsSequence])
        centroidFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).centroidList() for ts in obsSequence])


    centroidRa = '03:37:43.826'
    centroidDec = '14:15:14.828'
    haOffset = 150.
    print flatFilenames

    nWorkers = 7
    activeWorkers = []
    for iSeq,obsSeq in enumerate(obsFilenames):
        for iObs,obsPath in enumerate(obsSeq):
            print 'starting',obsPath
            centroidPath = centroidFilenames[iSeq][iObs]
#            xGuess = guessCols[iSeq]
#            yGuess = guessRows[iSeq]
            yGuess,xGuess = centerPixelGuesses[iSeq]
#            hotPath = timeMaskFilenames[iSeq][iObs]
#            flatPath = flatFilenames[iSeq]
            tstamp = obsSequences[iSeq][iObs]
            worker = multiprocessing.Process(target=centroidObs,args=(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess,savePath,tstamp))
            worker.start()
            activeWorkers.append(worker)
            if len(activeWorkers) >= nWorkers:
                for worker in activeWorkers:
                    worker.join()
                activeWorkers = []

            #centroidObs(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess,hotPath,flatPath)

    
