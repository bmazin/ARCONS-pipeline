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
from util.utils import confirm
from photometry.PSFphotometry import PSFphotometry


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

    nWorkers = 1
    activeWorkers = []
    print centerPixelGuesses
    bMakeNewHotPix = False
    integrationTime = 100
    exptime = 300
    for iSeq,obsSeq in enumerate(obsFilenames):
        for iObs,obsPath in enumerate(obsSeq):
            print 'starting',obsPath
            centroidPath = centroidFilenames[iSeq][iObs]
            yGuess,xGuess = centerPixelGuesses[iSeq]
            xyguess = [xGuess,yGuess]
            hotPath = timeMaskFilenames[iSeq][iObs]
            flatPath = flatFilenames[iSeq]

            obs = ObsFile(obsPath)
            print obsPath,obs.getFromHeader('exptime'),obs
            obs.loadBestWvlCalFile()
            obs.loadFlatCalFile(flatPath)
            obs.setWvlCutoffs(3000,8000)
            
            if bMakeNewHotPix:
                hp.findHotPixels(obsFile=obs,outputFileName=hotPath,display=True,fwhm=2.,boxSize=5, nSigmaHot=4.0,)
            obs.loadHotPixCalFile(hotPath,switchOnMask=False)
#switchOnHotPixTimeMask
#switchOffHotPixTimeMask
#cc.centroidCalc(obs,centroidRa,centroidDec,guessTime=300,integrationTime=100,secondMaxCountsForDisplay=2000,HA_offset=haOffset,xyapprox=[xGuess,yGuess],outputFileName=centroidPath,usePsfFit=True)

            for iFrame in range(0,exptime,integrationTime):

                beforeImgDict = obs.getPixelCountImage(firstSec=iFrame, integrationTime= integrationTime, weighted=True,fluxWeighted=False, getRawCount=False,scaleByEffInt=False)
                beforeImg=beforeImgDict['image']


                obs.switchOnHotPixTimeMask()
                afterImgDict = obs.getPixelCountImage(firstSec=iFrame, integrationTime= integrationTime, weighted=True,fluxWeighted=False, getRawCount=False,scaleByEffInt=False)
                afterImg=afterImgDict['image']

                afterImgDict2 = obs.getPixelCountImage(firstSec=iFrame, integrationTime= integrationTime, weighted=True,fluxWeighted=False, getRawCount=False,scaleByEffInt=True)
                afterImg2=afterImgDict2['image']
                afterImg2[np.invert(np.isfinite(afterImg2))] = 0.
                obs.switchOffHotPixTimeMask()

                plotArray(beforeImg,title='before')
                plotArray(afterImg,title='after')
                plotArray(afterImg2,title='after scaled')

                psfPhot = PSFphotometry(afterImg2,centroid=[xyguess],verbose=True)
                psfDict = psfPhot.PSFfit(aper_radius=6)
                fitModelImg = psfDict['fitModelImg']
                fitModelImg[afterImg2==0] = 0.
                errImg = np.sqrt(fitModelImg)
                diffImg = np.abs(fitModelImg-afterImg2)/errImg

                plotArray(diffImg)

                confirmResponse = confirm('Continue?')
                if  confirmResponse == False:
                    exit(1)

            del obs

    
