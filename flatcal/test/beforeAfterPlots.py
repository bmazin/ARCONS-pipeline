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
    flatCalDate = '20121211'
    flatCalTstamp = '20121212-074700'
    obsFN = FileName(run=run,date=date,tstamp=obsTimestamp)
    obsFileName = obsFN.obs()
    timeMaskFileName = obsFN.timeMask()
    wvlFileName = FileName(run=run,date=date,tstamp=wvlCal).calSoln()
    flatFileName = FileName(run=run,date=flatCalDate,tstamp=flatCalTstamp).flatSoln()
    
    if not os.path.exists(timeMaskFileName):
        print 'Running hotpix for ',obsFileName
        hp.findHotPixels(obsFileName,timeMaskFileName)
        print "Flux file pixel mask saved to %s"%(timeMaskFileName)

    obs = ObsFile(obsFileName)
    obs.loadTimeAdjustmentFile(FileName(run='PAL2012').timeAdjustments())
    obs.loadWvlCalFile(wvlFileName)
    obs.loadFlatCalFile(flatFileName)
    obs.loadHotPixCalFile(timeMaskFileName)

    #obs.setWvlCutoffs(4000,8000)

    #get image before and after flat cal
    print 'getting images'
    beforeImgDict = obs.getPixelCountImage(weighted=False,fluxWeighted=False,scaleByEffInt=True)

    rawCubeDict = obs.getSpectralCube(weighted=False)
    rawCube = np.array(rawCubeDict['cube'],dtype=np.double)
    effIntTime = rawCubeDict['effIntTime']
    maxIntTime = np.max(effIntTime)
    #add third dimension for broadcasting
    effIntTime3d = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
    rawCube *= maxIntTime / effIntTime3d
    rawCube[np.isnan(rawCube)] = 0
    rawCube[rawCube == np.inf] = 0
    beforeImg = np.sum(rawCube,axis=-1)
    print 'finished first cube'

    flatCubeDict = obs.getSpectralCube(weighted=True)
    flatCube = np.array(flatCubeDict['cube'],dtype=np.double)
    effIntTime = flatCubeDict['effIntTime']
    maxIntTime = np.max(effIntTime)
    #add third dimension for broadcasting
    effIntTime3d = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
    flatCube *= maxIntTime / effIntTime3d
    flatCube[np.isnan(flatCube)] = 0
    flatCube[flatCube == np.inf] = 0
    afterImg = np.sum(flatCube,axis=-1)

    plotArray(title='before flatcal',image=beforeImg)
    plotArray(title='after flatcal',image=afterImg)

    print 'before sdev',np.std(beforeImg[afterImg!=0])
    print 'after sdev',np.std(afterImg[afterImg!=0])

    np.savez('flatCubeGem.npz',beforeImg=beforeImg,afterImg=afterImg,rawCube=rawCube,flatCube=flatCube)
    

if __name__=='__main__':
    main()

