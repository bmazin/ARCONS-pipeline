#!/bin/python
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.popup import PopUp,plotArray
import util.utils
import matplotlib.pyplot as plt
import numpy as np
import datetime
import tables
import ephem
import matplotlib
import matplotlib.cm as cm
import os
import hotpix.hotPixels as hp

def main():

    #open the sky file for hr9087
    run = 'PAL2012'
    date = '20121210'
    wvlCal = '20121211-052230'
    obsTimestamp = '20121211-051650'
    flatCalDate = '20121210'
    obsFN = FileName(run=run,date=date,tstamp=obsTimestamp)
    obsFileName = obsFN.obs()
    timeMaskFileName = obsFN.timeMask()
    wvlFileName = FileName(run=run,date=date,tstamp=wvlCal).calSoln()
    flatFileName = FileName(run=run,date=flatCalDate,tstamp='').flatSoln()
    
    if not os.path.exists(timeMaskFileName):
        print 'Running hotpix for ',obsFileName
        hp.findHotPixels(obsFileName,timeMaskFileName)
        print "Flux file pixel mask saved to %s"%(timeMaskFileName)

    obs = ObsFile(obsFileName)
    obs.loadTimeAdjustmentFile(FileName(run='PAL2012').timeAdjustments())
    obs.loadWvlCalFile(wvlFileName)
    obs.loadFlatCalFile(flatFileName)
    obs.loadHotPixCalFile(timeMaskFileName)

    obs.setWvlCutoffs(4000,8000)

    #get image before and after flat cal
    print 'getting images'
    beforeImgDict = obs.getPixelCountImage(weighted=False,fluxWeighted=False,scaleByEffInt=True)
    beforeImg = beforeImgDict['image']
    beforeImg[np.isnan(beforeImg)]=0
    beforeImg[beforeImg == np.inf] = 0
    plotArray(title='before flatcal',image=beforeImg)
    afterImg = obs.getPixelCountImage(weighted=True,fluxWeighted=False,scaleByEffInt=True)['image']
    afterImg[np.isnan(afterImg)]=0
    afterImg[afterImg == np.inf] = 0
    print 'plotting images'
    plotArray(title='after flatcal',image=afterImg)

    print 'before sdev',np.std(beforeImg[afterImg!=0])
    print 'after sdev',np.std(afterImg[afterImg!=0])

    np.savez('beforeAfterImgs.npz',beforeImg=beforeImg,afterImg=afterImg)
    

if __name__=='__main__':
    main()

