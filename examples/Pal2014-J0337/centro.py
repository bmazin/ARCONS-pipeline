import numpy as np
import tables
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
import os
import hotpix.hotPixels as hp
from util.FileName import FileName
from util.popup import plotArray
import astrometry.CentroidCalc as cc
from photonlist.PhotList import writePhotonList

if __name__=='__main__':
    run = 'PAL2014'
    date = '20140924'
    tstamp = '20140925-112538'

    obsFN = FileName(run=run,date=date,tstamp=tstamp)
    obsPath = obsFN.obs()
     
    flatPath = FileName(run=run,date=date).flatSoln()
    hotPath = obsFN.timeMask()
    centroidPath = obsFN.centroidList()
    
    obs = ObsFile(obsPath)
    if not os.path.exists(hotPath):
        hp.findHotPixels(obsPath,hotPath)
    obs.loadHotPixCalFile(hotPath,switchOnMask=False)
    obs.loadBestWvlCalFile()
    obs.loadFlatCalFile(flatPath)
    obs.setWvlCutoffs(3000,8000)

    
    centroidRa = '03:37:43.826'
    centroidDec = '14:15:14.828'
    haOffset = 150.

    imgDict = obs.getPixelCountImage(integrationTime=60,scaleByEffInt=True)
    plotArray(imgDict['image'])
    

    xGuess = 18 #col
    yGuess = 15 #row
    cc.centroidCalc(obs,centroidRa,centroidDec,guessTime=300,integrationTime=30,secondMaxCountsForDisplay=2000,HA_offset=haOffset,xyapprox=[xGuess,yGuess],outputFileName=centroidPath)

    
