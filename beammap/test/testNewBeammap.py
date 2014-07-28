'''
Author: Matt Strader                    Date: May 21, 2014
Test if a remapped beammap is actually remapped
'''
from util.ObsFile import ObsFile
from util.FileName import FileName
from util.popup import PopUp,plotArray
import util.utils
import matplotlib.pyplot as plt
import numpy as np
import tables
import os
from beammap import remapPixels
import unittest

def testLoadBeammap():
    '''
    Test if a remapped beammap is actually remapped
    '''
    #open an obs file from PAL2012,the sky file for hr9087
    #we'll use the embedded beammap file, which has some misplaced pixels
    run = 'PAL2012'
    date = '20121210'
    obsTimestamp = '20121211-051650'
    obsFN = FileName(run=run,date=date,tstamp=obsTimestamp)
    obsFileName = obsFN.obs()
    obs = ObsFile(obsFileName)
    
    beammapFileName = obsFN.beammap()
    
    #load the PixelMap for PAL2012 to know which pixels should be remapped
    pixMap = remapPixels.PixelMap(obsFN.pixRemap())
    pixMapSourceList,pixMapDestList = pixMap.getRemappedPix()
    
    #load the corrected beammap into the obs
    obs.loadBeammapFile(beammapFileName)

    #check that each pixel that should be moved is moved
    #by comparing the embedded beammap and the loaded corrected one
    for source,dest in zip(pixMapSourceList,pixMapDestList):
        assert obs.beamImage[dest] == obs.file.root.beammap.beamimage[source]

    obs.file.close()

if __name__ == '__main__':
    testLoadBeammap()
