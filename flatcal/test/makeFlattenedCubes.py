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
from astropy.stats.funcs import sigma_clip

if __name__=='__main__':

    np.seterr(all='ignore')
    outFolder = '/Scratch/dataProcessing/flatTests/'
    #open the sky file for hr9087
    run = 'PAL2012'
    wvlCal = '' #'20121211-052230'

    flatCalDate = '20121211'
    #flatCalDate = '20140924'
    #flatCalTstamp = '20140925dome'
    #flatCalTstamp = '20121212-074700'
    flatCalTstamp = ''
    flatOutputLabel = flatCalDate
    lowerWvlCut = 4000
    upperWvlCut = 8000

    #date = '20121210'
    #obsTimestamp = '20121211-051650' #hr9087
    #obsTimestamp = '20121211-124809' #J0926
    #obsTimestamp = '20121211-125312' #J0926
    #obsTimestamp = '20121211-125814' #J0926
    #obsTimestamp = '20121211-130316' #J0926
    #obsTimestamp = '20121211-130818' #J0926
    #obsTimestamp = '20121211-131320' #J0926
    #obsTimestamp = '20121211-131822' #J0926
    #obsTimestamp = '20121211-132324' #last J0926

    #obsTimestamp = '20121211-134751' # feige 66 sky
    #obsTimestamp = '20121211-134812' # feige 66 sky
    #obsTimestamp = '20121211-134914' # feige 66 sky
    #obsTimestamp = '20121211-135016' #early twilight
    #obsTimestamp = '20121211-135118' #later twilight
    #obsTimestamp = '20121211-135220' #later twilight
    #obsTimestamp = '20121211-135322' #later twilight
    #obsTimestamp = '20121211-135424' #later twilight
    #obsTimestamp = '20121211-135526' #later twilight
    #obsTimestamp = '20121211-135628' #late twilight
    #obsTimestamp = '20121211-135730' #late twilight


    date = '20121211'
    #obsTimestamp = '20121212-085730' #geminga
    #obsTimestamp = '20121212-090233' #geminga
    #obsTimestamp = '20121212-090735' #geminga
    #obsTimestamp = '20121212-091237' #geminga
    
    #obsTimestamp = '20121212-095809' #ptfo 8-8695
    #obsTimestamp = '20121212-100311' #ptfo 8-8695
    #obsTimestamp = '20121212-100813' #ptfo 8-8695
    #obsTimestamp = '20121212-101315' #ptfo 8-8695
    #obsTimestamp = '20121212-101817' #ptfo 8-8695

    #obsTimestamp = '20121212-103030' #geminga
    #obsTimestamp = '20121212-104035' #geminga
    #obsTimestamp = '20121212-105326' #geminga
    #obsTimestamp = '20121212-110330' #geminga
    obsTimestamp = '20121212-111334' #geminga

    #obsTimestamp = '20121212-112709' #J0926
    #obsTimestamp = '20121212-113212' #J0926
    #obsTimestamp = '20121212-113714' #J0926
    #obsTimestamp = '20121212-114216' #J0926
    ##obsTimestamp = '20121212-114718' #J0926
    #obsTimestamp = '20121212-115220' #J0926
    #obsTimestamp = '20121212-115722' #J0926
    #obsTimestamp = '20121212-120224' #J0926
    #obsTimestamp = '20121212-120727' #J0926
    #obsTimestamp = '20121212-121229' #J0926
    #obsTimestamp = '20121212-121732' #J0926
    #obsTimestamp = '20121212-122234' #J0926
    #obsTimestamp = '20121212-122736' #J0926
    #obsTimestamp = '20121212-123238' #J0926
    #obsTimestamp = '20121212-123740' #J0926
    #obsTimestamp = '20121212-124242' #J0926
    #obsTimestamp = '20121212-124744' #J0926
    #obsTimestamp = '20121212-125246' #J0926
    #obsTimestamp = '20121212-125748' #J0926
    #obsTimestamp = '20121212-130250' #J0926
    #obsTimestamp = '20121212-130752' #J0926
    #obsTimestamp = '20121212-131254' #J0926
    #obsTimestamp = '20121212-131756' #J0926
    #obsTimestamp = '20121212-132258' #J0926
    #obsTimestamp = '20121212-132800' #J0926
    #obsTimestamp = '20121212-133303' #J0926
    #obsTimestamp = '20121212-132258' #J0926
    #obsTimestamp = '20121212-132800' #J0926
    #obsTimestamp = '20121212-133303' #J0926
    ##obsTimestamp = '20121212-134024' #early twilight
    #obsTimestamp = '20121212-134127' #twilight
    #obsTimestamp = '20121212-134229' #twilight
    #obsTimestamp = '20121212-134331' #twilight
    #obsTimestamp = '20121212-134433' #twilight
    ##obsTimestamp = '20121212-134535' #twilight
    #obsTimestamp = '20121212-134637' #twilight
    #obsTimestamp = '20121212-134739' #twilight
    #obsTimestamp = '20121212-134841' #twilight
    #obsTimestamp = '20121212-134943' #twilight
    #obsTimestamp = '20121212-135045' #twilight
    #obsTimestamp = '20121212-135147' #twilight
    #obsTimestamp = '20121212-135249' #twilight
    #obsTimestamp = '20121212-135351' #twilight
    #obsTimestamp = '20121212-135453' #twilight
    #obsTimestamp = '20121212-135555' #twilight
    #obsTimestamp = '20121212-135657' #twilight
    #obsTimestamp = '20121212-135759' #twilight
    #obsTimestamp = '20121212-135901' #twilight
    ##obsTimestamp = '20121212-140003' #twilight
    ##obsTimestamp = '20121212-140105' #late twilight


    #date = '20140923'
    #obsTimestamp = '20140924-065535' #early 1SWASP sky
    #obsTimestamp = '20140924-122405' #last 1SWASP
    #obsTimestamp = '20140924-123738' #early twilight


    #date = '20140924'
    #obsTimestamp ='20140925-030054' #G24-9 sky, faint objects on array
    #obsTimestamp ='20140925-052248' #V407 Vul
    #obsTimestamp ='20140925-052804' #V407 Vul
    #obsTimestamp ='20140925-053307' #V407 Vul
    #obsTimestamp ='20140925-053810' #V407 Vul
    #obsTimestamp ='20140925-110513' #psr J0337
    #obsTimestamp ='20140925-111029' #psr J0337
    #obsTimestamp ='20140925-111532' #psr J0337
    #obsTimestamp ='20140925-112035' #psr J0337
    #obsTimestamp ='20140925-112538' #psr J0337
    #obsTimestamp ='20140925-113041' #psr J0337
    #obsTimestamp ='20140925-113544' #psr J0337
    #obsTimestamp ='20140925-114047' #psr J0337
    #obsTimestamp ='20140925-114550' #psr J0337
    #obsTimestamp ='20140925-115053' #mid psr J0337
    #obsTimestamp ='20140925-115605' #psr J0337
    #obsTimestamp ='20140925-120112' #psr J0337
    #obsTimestamp ='20140925-120615' #psr J0337
    #obsTimestamp ='20140925-121118' #psr J0337
    #obsTimestamp ='20140925-121621' #psr J0337
    #obsTimestamp ='20140925-122124' #psr J0337
    #obsTimestamp ='20140925-122627' #psr J0337
    #obsTimestamp ='20140925-123130' #psr J0337
    #obsTimestamp ='20140925-123633' #last PSR J0337 file for 20140924
    #obsTimestamp ='20140925-124254' #early twilight
    #obsTimestamp ='20140925-124357' #early twilight
    #obsTimestamp ='20140925-124500' #early twilight
    #obsTimestamp ='20140925-124603' #early twilight
    #obsTimestamp ='20140925-124706' #early twilight
    #obsTimestamp ='20140925-124809' #early twilight
    #obsTimestamp ='20140925-124912' #early twilight
    #obsTimestamp ='20140925-125015' #late twilight
    #obsTimestamp ='20140925-125118' #late twilight
    #obsTimestamp ='20140925-125221' #late twilight
    #obsTimestamp ='20140925-125324' #late twilight
    #obsTimestamp ='20140925-125427' #late twilight

    #date = '20140925'
    #obsTimestamp ='20140926-104528' #psr J0337


    obsFN = FileName(run=run,date=date,tstamp=obsTimestamp)
    obsFileName = obsFN.obs()
    print obsFileName
    timeMaskFileName = obsFN.timeMask()
    if wvlCal != '':
        wvlFileName = FileName(run=run,date=date,tstamp=wvlCal).calSoln()
    flatFileName = FileName(run=run,date=flatCalDate,tstamp=flatCalTstamp).flatSoln()
    illumFileName = FileName(run=run,date=flatCalDate,tstamp=flatCalTstamp).illumSoln()
    
    if not os.path.exists(timeMaskFileName):
        print 'Running hotpix for ',obsFileName
        hp.findHotPixels(inputFileName=obsFileName,outputFileName=timeMaskFileName)
        print "Flux file pixel mask saved to %s"%(timeMaskFileName)


    obs = ObsFile(obsFileName)
    exptime = obs.getFromHeader('exptime')
    #obs.loadTimeAdjustmentFile(FileName(run='PAL2012').timeAdjustments())
    if wvlCal != '':
        obs.loadWvlCalFile(wvlFileName)
    else:
        obs.loadBestWvlCalFile()
    obs.loadHotPixCalFile(timeMaskFileName)

    obs.setWvlCutoffs(lowerWvlCut,upperWvlCut)

    obs.loadFlatCalFile(flatFileName)
    #get image before and after flat cal
    #print 'getting images'
    rawCubeDict = obs.getSpectralCube(weighted=False)
    rawCube = np.array(rawCubeDict['cube'],dtype=np.double)
    effIntTime = rawCubeDict['effIntTime']
    maxIntTime = np.max(effIntTime)
    #add third dimension for broadcasting
    effIntTime3d = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
    rawCube *= maxIntTime / effIntTime3d
    rawImg = np.sum(rawCube,axis=-1)
    rawImg[rawImg == 0] = np.nan
    rawImg[~np.isfinite(rawImg)] = np.nan
    #print 'finished raw cube'

    flatCubeDict = obs.getSpectralCube(weighted=True)
    flatCube = np.array(flatCubeDict['cube'],dtype=np.double)
    effIntTime = flatCubeDict['effIntTime']
    maxIntTime = np.max(effIntTime)
    #add third dimension for broadcasting
    effIntTime3d = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
    flatCube *= maxIntTime / effIntTime3d
    flatImg = np.sum(flatCube,axis=-1)
    flatImg[flatImg == 0] = np.nan
    flatImg[~np.isfinite(flatImg)] = np.nan
    #print 'finished flat cube'

    #Now load the illum cal and get another flattened cube
    obs.loadFlatCalFile(illumFileName)
    illumCubeDict = obs.getSpectralCube(weighted=True)
    illumCube = np.array(illumCubeDict['cube'],dtype=np.double)
    effIntTime = illumCubeDict['effIntTime']
    maxIntTime = np.max(effIntTime)
    #add third dimension for broadcasting
    effIntTime3d = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
    illumCube *= maxIntTime / effIntTime3d
    illumImg = np.sum(illumCube,axis=-1)
    illumImg[illumImg == 0] = np.nan
    illumImg[~np.isfinite(illumImg)] = np.nan
    #print 'finished illum cube'

    #plotArray(title='raw',image=rawImg)
    #plotArray(title='with flatcal',image=flatImg)
    #plotArray(title='with illumcal',image=illumImg)

    #print 'raw sdev',np.std(rawImg[rawImg!=0])
    #print 'flat sdev',np.std(flatImg[flatImg!=0])
    #print 'illum sdev',np.std(illumImg[illumImg!=0])

    np.savez(os.path.join(outFolder,'flattenedCubes_obs{}_flat{}_wvl{}-{}.npz'.format(obsTimestamp,flatOutputLabel,lowerWvlCut,upperWvlCut)),rawImg=rawImg,flatImg=flatImg,illumImg=illumImg,rawCube=rawCube,flatCube=flatCube,illumCube=illumCube,exptime=exptime)
    


