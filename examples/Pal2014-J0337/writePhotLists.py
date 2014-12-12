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

def writePL(tstamp,obsPath,hotPath,flatPath,fluxPath):
    obs = ObsFile(obsPath)
    if not os.path.exists(hotPath):
        hp.findHotPixels(obsFile=obs,outputFileName=hotPath)
    obs.loadHotPixCalFile(hotPath,switchOnMask=False)
    obs.loadBestWvlCalFile()
    obs.loadFlatCalFile(flatPath)
    obs.loadFluxCalFile(fluxPath)
    obs.setWvlCutoffs(3000,11000)
    obs.loadCentroidListFile(centroidPath)
    writePhotonList(obs,photListDescription=PulsarPhotonList)
    del obs

if __name__=='__main__':
    savePath = '/Scratch/dataProcessing/J0337/'
    run = 'PAL2014'
    sunsetDates = []
    flatDates = []
    obsSequences = []
    parFiles = []
    

    obsSequence0 = [
    '20140925-111532',
    '20140925-112035',
    '20140925-112538',
    '20140925-113041',
    '20140925-113544',
    '20140925-114047',
    '20140925-114550',
    '20140925-115053',
    '20140925-115605',
    '20140925-120112',
    '20140925-120615',
    '20140925-121118',
    '20140925-121621',
    '20140925-122124',
    '20140925-122627',
    '20140925-123130',
    '20140925-123633'
    ]
    sunsetDates.append('20140924')
    flatDates.append('20140924')
    obsSequences.append(obsSequence0)
    parFiles.append('J0337-mjd-56924.par')
    

    obsSequence1 = [
    '20140926-081832',
    '20140926-082336',
    '20140926-082840',
    '20140926-083344',
    '20140926-083848',
    '20140926-084352',
    '20140926-084856',
    '20140926-085400',
    '20140926-085904',
    '20140926-090408',
    '20140926-090912',
    '20140926-091416',
    '20140926-091920',
    '20140926-092424',
    '20140926-092928',
    '20140926-093432',
    '20140926-093936',
    '20140926-094440',
    '20140926-094944',
    '20140926-095448',
    '20140926-095952',
    '20140926-100456',
    '20140926-101000',
    '20140926-101504',
    '20140926-102008',
    '20140926-102512',
    '20140926-103016',
    '20140926-103520',
    '20140926-104024',
    '20140926-104528'
    ]
    sunsetDates.append('20140925')
    flatDates.append('20140924')
    obsSequences.append(obsSequence1)
    parFiles.append('J0337-mjd-56924.par')

    obsSequence2 = [
    '20141021-063246',
    '20141021-063750',
    '20141021-064253',
    '20141021-064758',
    '20141021-065303',
    '20141021-065807',
    '20141021-070311',
    '20141021-071138',
    '20141021-071656',
    '20141021-072200',
    '20141021-072704',
    '20141021-073208',
    '20141021-073713',
    '20141021-074216',
    '20141021-074720',
    '20141021-075225',
    '20141021-075729',
    '20141021-080232',
    '20141021-080737',
    '20141021-081242',
    '20141021-081747',
    '20141021-082252',
    '20141021-082755',
    '20141021-083300',
    '20141021-083805',
    '20141021-084310',
    '20141021-084814',
    '20141021-085318',
    '20141021-085822',
    '20141021-090326',
    '20141021-090832',
    '20141021-091336',
    '20141021-091840',
    '20141021-092345',
    '20141021-092851'
    ]
    sunsetDates.append('20141020')
    flatDates.append('20141020')
    obsSequences.append(obsSequence2)
    parFiles.append('J0337-mjd-56951.par')


    obsFilenames = []
    timeMaskFilenames = []
    flatFilenames = []
    centroidFilenames = []
    plPaths = []
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        sunsetDate = sunsetDates[iSeq]
        obsFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).obs() for ts in obsSequence])
        timeMaskFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).timeMask() for ts in obsSequence])
        #fluxFilenames.append(FileName(run=run,date=fluxCalDates[iSeq],tstamp=fluxCals[iSeq]).fluxSoln())
        flatFilenames.append(FileName(run=run,date=flatDates[iSeq],tstamp='').flatSoln())
        #plFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).crabList() for ts in obsSequence])
        centroidFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).centroidList() for ts in obsSequence])
        plPaths.append([FileName(run=run,date=sunsetDate,tstamp=ts).photonList() for ts in obsSequence])

    fluxPath = FileName(run='PAL2012',date='20121211',tstamp='absolute_021727').fluxSoln()


    nWorkers = 7
    activeWorkers = []
    for iSeq,obsSeq in enumerate(obsFilenames):
        for iObs,obsPath in enumerate(obsSeq):
            print 'starting',obsPath
            centroidPath = centroidFilenames[iSeq][iObs]
            hotPath = timeMaskFilenames[iSeq][iObs]
            flatPath = flatFilenames[iSeq]

            tstamp = obsSequences[iSeq][iObs]
            worker = multiprocessing.Process(target=writePL,args=(tstamp,obsPath,hotPath,flatPath,fluxPath))
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

    workingDir = '/Scratch/dataProcessing/J0337/'
    for iSeq,plSeq in enumerate(plPaths):
        parFile = parFiles[iSeq]
        for iPL,plPath in enumerate(plSeq):
            photonTiming.timePhotonList(plPath,parFile=parFile,bPulsarTiming=True,timingProgram='tempo',verbose=True,nPhotonsPerProcess=1e5,workingDir=workingDir)
    
    print 'done'
            #centroidObs(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess,hotPath,flatPath)

    
