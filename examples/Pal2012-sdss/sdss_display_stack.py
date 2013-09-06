#!/bin/python

### see 0926params.dict to change variables. This program creates an image stack of the data which is then input into sdssfitpsf.py. There a 2-d gaussian is fit to each image, and the resulting file is used to create a lightcurve using sdsslightcurve.py, or plotFullData.py

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
param.read_from_file('0926params.dict')

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

wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]

#wvlCalFilenames = ['/Scratch/waveCalSolnFiles/%s/calsol_%s.h5'%(sunsetDate,sunsetDate) for sunsetDate in sunsetDates]

flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]

#wvlCalFilenames = ['/Scratch/waveCalSolnFiles/%s/wvlCalFiles_combined_Dec%s2012.h5'%(sunsetDate,sunsetDate[6:]) for sunsetDate in sunsetDates]


#flatCalFilenames = ['/Scratch/flatCalSolnFiles/%s/flatsol_%s.h5'%(sunsetDate,sunsetDate) for sunsetDate in sunsetDates]

#use for Dec8
flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'
flatCalFilenames[1] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

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
        obsFn = FileName(run=run,date=sunsetDate,tstamp=ts).obs()
        ob = ObsFile(obsFn)
        ob.loadTimeAdjustmentFile(FileName(run=run).timeAdjustments())
        index1 = obsFn.find('_')
        hotPixFn = '/Scratch/timeMasks/timeMask' + obsFn[index1:]
        if not os.path.exists(hotPixFn):
            hp.findHotPixels(obsFn,hotPixFn)
            print "Flux file pixel mask saved to %s"%(hotPixFn)
        ob.loadHotPixCalFile(hotPixFn,switchOnMask=True)
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
        ob.setWvlCutoffs(wvlLowerCutoff,wvlUpperCutoff)

        bad_solution_mask=np.zeros((46,44))
        bad_count=0;
        for y in range(46):
            for x in range(44):
                if (ob.wvlRangeTable[y][x][1] < wvlUpperCutoff) or (ob.wvlRangeTable[y][x][0] > wvlLowerCutoff):
                    bad_solution_mask[y][x] = 1

        unix = ob.getFromHeader('unixtime')
        startJD = unix/86400.+2440587.5
        nSecInFile = ob.getFromHeader('exptime')
        #tic = time()
        deadMask = ob.getDeadPixels()

	#print 'Dead mask load time = ', time()-tic
        for sec in np.arange(0,nSecInFile,integrationTime):
            jd = startJD + sec/(24.*3600.) + integrationTime/2./(24.*3600.)#add seconds offset to julian date, move jd to center of bin
            print count,jd
            count+=1
            times.append(jd)
            titles.append('%.6f'%jd)
            frameData = ob.getPixelCountImage(firstSec=sec,integrationTime=integrationTime,weighted=True,scaleByEffInt=True)
            frame = frameData['image']         
            showFrame = np.array(frame)
            showframes.append(showFrame)
            frame[deadMask == 0] = np.nan
            #frame[bad_solution_mask == 1] = np.nan
            frames.append(frame)
          
cube = np.dstack(frames)
times = np.array(times)

np.savez(npzFileName,stack=cube,jd=times)
print 'saved'
utils.makeMovie(showframes,cbar=True,frameTitles=titles,outName=gifFileName)

