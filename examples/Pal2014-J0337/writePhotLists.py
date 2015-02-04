import numpy as np
import tables
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
import os
import hotpix.hotPixels as hp
from util.FileName import FileName
from util.popup import plotArray,PopUp
import astrometry.CentroidCalc as cc
import multiprocessing
from photonlist.photlist import writePhotonList
from headers.ArconsHeaders import PulsarPhotonList
from timing import photonTiming
from util.readDict import readDict

def writePL(tstamp,obsPath,centroidPath):
    obs = ObsFile(obsPath)
    obs.loadAllCals()
    obs.setWvlCutoffs(3000,11000)
    obs.loadCentroidListFile(centroidPath)
    writePhotonList(obs,photListDescription=PulsarPhotonList,checkForExisting=False)
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
    plPaths = []
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        sunsetDate = sunsetDates[iSeq]
        obsFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).obs() for ts in obsSequence])
        centroidFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).centroidList() for ts in obsSequence])
        plPaths.append([FileName(run=run,date=sunsetDate,tstamp=ts).photonList() for ts in obsSequence])


    nWorkers = 7
    activeWorkers = []
    for iSeq,obsSeq in enumerate(obsFilenames):
        for iObs,obsPath in enumerate(obsSeq):
            print 'starting',obsPath
            centroidPath = centroidFilenames[iSeq][iObs]

            tstamp = obsSequences[iSeq][iObs]
            worker = multiprocessing.Process(target=writePL,args=(tstamp,obsPath,centroidPath))
            worker.start()
            activeWorkers.append(worker)
            if len(activeWorkers) >= nWorkers:
                for worker in activeWorkers:
                    worker.join()
                activeWorkers = []

    if len(activeWorkers) > 0:
        for worker in activeWorkers:
            worker.join()
        activeWorkers = []

#    workingDir = '/Scratch/dataProcessing/J0337/'
#    for iSeq,plSeq in enumerate(plPaths):
#        parFile = parFiles[iSeq]
#        for iPL,plPath in enumerate(plSeq):
#            photonTiming.timePhotonList(plPath,parFile=parFile,bPulsarTiming=True,timingProgram='tempo',verbose=True,nPhotonsPerProcess=1e5,workingDir=workingDir)
    
    print 'done'

    
