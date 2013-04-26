#!/bin/python

'''Julian's version of sdss_display_stack.py - for now, just for trying to get some
photometry out on CoRoT 18. In future will probably modify/rewrite as a generalised
photometry code'''

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util.readDict import readDict
from util import utils
import tables
import matplotlib.pyplot as plt
import hotpix.hotPixels as hp
import os
from time import time

param = readDict()
param.read_from_file('photometryParams.dict')

mkidDataDir = param['mkidDataDir']
intermDir = param['intermDir']

run = param['run']
flatSunsetLocalDate = param['flatSunsetLocalDate']
flatTimestamp = param['flatTimestamp']

wvlCalSunsetLocalDate = param['wvlCalSunsetLocalDate']
wvlCalTimestamp = param['wvlCalTimestamp']

expTime = param['expTime']
integrationTime = param['integrationTime']

utcDates = param['utcDates']
sunsetDates = param['sunsetDates']

calTimestamps = param['calTimestamps']
seqs = param['seqs']

wvlLowerCutoff = param['wvlLowerCutoff']
wvlUpperCutoff = param['wvlUpperCutoff']
npzFileName = param['npzFileName']
gifFileName = param['gifFileName']

timestampLists = [[utcDate+'-'+str(ts) for ts in seq] for utcDate,seq in zip(utcDates,seqs)]

wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp,mkidDataDir=mkidDataDir,
                            intermDir=intermDir).calSoln() 
                   for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]

##wvlCalFilenames[0] = '/Scratch/waveCalSolnFiles/20121210/calsol_20121211-074031.h5'
##wvlCalFilenames[1] = '/home/danica/optimusP/testing/forMatt/calsol_20121211-044853.h5'
##flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(['20121210','20121210'],['20121211-074031','20121211-074031'])]

flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp,mkidDataDir=mkidDataDir,
                             intermDir=intermDir).flatSoln()
                     for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]

#use for Dec8
#flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'
#flatCalFilenames[1] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

frames = []
showframes = []
times = []
titles = []
count= 0

NumFiles = []
for seq in seqs:
    NumFiles.append(len(seq))
NumFiles = sum(NumFiles)
print (NumFiles)*expTime/integrationTime,'frames to make'

for iSeq in range(len(seqs)):
    timestampList = timestampLists[iSeq]
    print timestampList
    wfn = wvlCalFilenames[iSeq]
    ffn = flatCalFilenames[iSeq]
    sunsetDate = sunsetDates[iSeq]
    
    for i,ts in enumerate(timestampList):
        print 'loading',ts
        obsFn = FileName(run=run,date=sunsetDate,tstamp=ts,mkidDataDir=mkidDataDir,
                         intermDir=intermDir).obs()
        ob = ObsFile(obsFn)
        ob.loadTimeAdjustmentFile(FileName(run=run,mkidDataDir=mkidDataDir,
                            intermDir=intermDir).timeAdjustments())
        
        #index1 = obsFn.find('_')
        #hotPixFn = '/Scratch/timeMasks/timeMask' + obsFn[index1:]
        hotPixFn = FileName(run=run,date=sunsetDate,tstamp=ts,mkidDataDir=mkidDataDir,
                         intermDir=intermDir).timeMask()
        
        if not os.path.exists(hotPixFn):
            hp.findHotPixels(obsFn,hotPixFn)
            print "Flux file pixel mask saved to %s"%(hotPixFn)
        
        print 'Loading cal files'
        ob.loadHotPixCalFile(hotPixFn,switchOnMask=True)
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
        ob.setWvlCutoffs(wvlLowerCutoff,wvlUpperCutoff)

        print 'Finding bad wvl pixels'
        bad_solution_mask=np.zeros((46,44))
        bad_count=0;
        for y in range(46):
            for x in range(44):
                if (ob.wvlRangeTable[y][x][1] < 11000):
                    bad_solution_mask[y][x] = 1
                    bad_count = bad_count+1

        startJD = ob.getFromHeader('jd')
        nSecInFile = ob.getFromHeader('exptime')
        #tic = time()
        print 'Finding dead pixels'
        deadMask = ob.getDeadPixels()
        #print 'Dead mask load time = ', time()-tic
        print 'Getting images'
        for sec in range(0,nSecInFile,integrationTime):
            jd = startJD + sec/(24.*3600.)#add seconds offset to julian date
            print count,jd
            count+=1
            times.append(jd)
            titles.append('%.6f'%jd)
            frameData = ob.getPixelCountImage(firstSec=sec,integrationTime=integrationTime,weighted=True)
            frame = frameData['image']         
            showFrame = np.array(frame)
            showframes.append(showFrame)
            frame[deadMask == 0] = np.nan
            frames.append(frame)
          
cube = np.dstack(frames)
times = np.array(times)

np.savez(npzFileName,stack=cube,jd=times)
print 'saved'
utils.makeMovie(showframes,cbar=True,frameTitles=titles,outName=gifFileName)

