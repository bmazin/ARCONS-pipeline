import numpy as np
import tables
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
import os
import hotpix.hotPixels as hp
from util.FileName import FileName
from util.popup import plotArray
import astrometry.CentroidCalc as cc
import multiprocessing

def centroidObs(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess,hotPath,flatPath):
    obs = ObsFile(obsPath)
    print obsPath,obs.getFromHeader('exptime'),obs
    if not os.path.exists(hotPath):
        hp.findHotPixels(obsFile=obs,outputFileName=hotPath)
    obs.loadHotPixCalFile(hotPath,switchOnMask=False)
    obs.loadBestWvlCalFile()
    obs.loadFlatCalFile(flatPath)
    obs.setWvlCutoffs(3000,8000)
    cc.centroidCalc(obs,centroidRa,centroidDec,guessTime=300,integrationTime=30,secondMaxCountsForDisplay=2000,HA_offset=haOffset,xyapprox=[xGuess,yGuess],outputFileName=centroidPath)
    print 'done centroid',centroidPath
    del obs

if __name__=='__main__':
    run = 'PAL2014'
    sunsetDates = []
    flatDates = []
    obsSequences = []
    guessRows = []
    guessCols = []
    

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
#    sunsetDates.append('20140924')
#    flatDates.append('20140924')
#    obsSequences.append(obsSequence0)
#    guessRows.append(16)
#    guessCols.append(18)
    

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
    '20140926-104528',
    '20140926-105032'
    ]
    sunsetDates.append('20140925')
    flatDates.append('20140924')
    obsSequences.append(obsSequence1)
    guessRows.append(16)
    guessCols.append(12)

    obsSequence2 = [
    '20141021-063246',
    '20141021-063750',
    '20141021-064253',
    '20141021-064758',
    '20141021-065303',
    '20141021-065807',
    '20141021-070311',
    '20141021-070814',
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
    '20141021-092851',
    '20141021-093356'
    ]
    sunsetDates.append('20141020')
    flatDates.append('20141020')
    obsSequences.append(obsSequence2)
    guessRows.append(11)
    guessCols.append(18)


    obsFilenames = []
    timeMaskFilenames = []
    flatFilenames = []
    centroidFilenames = []
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        sunsetDate = sunsetDates[iSeq]
        obsFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).obs() for ts in obsSequence])
        timeMaskFilenames.append([FileName(run=run,date=sunsetDate,tstamp=ts).timeMask() for ts in obsSequence])
        #fluxFilenames.append(FileName(run=run,date=fluxCalDates[iSeq],tstamp=fluxCals[iSeq]).fluxSoln())
        flatFilenames.append(FileName(run=run,date=flatDates[iSeq],tstamp='').flatSoln())
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
            xGuess = guessCols[iSeq]
            yGuess = guessRows[iSeq]
            hotPath = timeMaskFilenames[iSeq][iObs]
            flatPath = flatFilenames[iSeq]
            worker = multiprocessing.Process(target=centroidObs,args=(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess,hotPath,flatPath))
            worker.start()
            activeWorkers.append(worker)
            if len(activeWorkers) >= nWorkers:
                for worker in activeWorkers:
                    worker.join()
                activeWorkers = []

            #centroidObs(obsPath,centroidPath,centroidRa,centroidDec,haOffset,xGuess,yGuess,hotPath,flatPath)

    
