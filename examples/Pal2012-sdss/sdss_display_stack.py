#!/bin/python

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util.readDict import readDict
from util import utils
import tables
import matplotlib.pyplot as plt
from hotpix import hotPixels

params = readDict()
params.read_from_file('0926params.dict')

run = params['run']
flatSunsetLocalDate = params['flatSunsetLocalDate']
flatTimestamp = params['flatTimestamp']

wvlCalSunsetLocalDate = params['wvlCalSunsetLocalDate']
wvlCalTimestamp = params['wvlCalTimestamp']

expTime = params['expTime']
integrationTime = params['integrationTime']

utcDates = params['utcDates']
sunsetDates = params['sunsetDates']

calTimestamps = params['calTimestamps']
seqs = params['seqs']

wvlLowerCutoff = params['wvlLowerCutoff']
wvlUpperCutoff = params['wvlUpperCutoff']
npzFileName = params['npzFileName']
gifFileName = params['gifFileName']

timestampLists = [[utcDate+'-'+str(ts) for ts in seq] for utcDate,seq in zip(utcDates,seqs)]

wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
#wvlCalFilenames[0] = '/Scratch/waveCalSolnFiles/20121210/calsol_20121211-074031.h5'
#wvlCalFilenames[1] = '/home/danica/optimusP/testing/forMatt/calsol_20121211-044853.h5'
#flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(['20121210','20121210'],['20121211-074031','20121211-074031'])]

flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
#flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'
#flatCalFilenames[1] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

#/Scratch/waveCalSolnFiles/20121208/calsol_20121209-131132.h5

frames = []
showframes = []
times = []
titles = []
#plt.ion()
count= 0

NumFiles = []
for seq in seqs:
    NumFiles.append(len(seq))
NumFiles = sum(NumFiles)

print (NumFiles)*expTime/integrationTime,'frames to make'

#print (len(seqs))*expTime/integrationTime,'frames to make'
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
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
        ob.setWvlCutoffs(wvlLowerCutoff,wvlUpperCutoff)

	bad_solution_mask=np.zeros((46,44))
	bad_count=0;
	for y in range(46):
	    for x in range(44):
		if (ob.wvlRangeTable[y][x][1] < 11000):
		    bad_solution_mask[y][x] = 1
		    bad_count = bad_count+1
        print bad_count

#        row1 = 19
#        col1 = 30
#        row2 = 8
#        col2 = 30
#        print '%d,%d'%(row1,col1),ob.getPixelCount(row1,col1,weighted=False,integrationTime=30),ob.getPixelCount(row1,col1,weighted=True,integrationTime=30)
#        print '%d,%d'%(row2,col2),ob.getPixelCount(row2,col2,weighted=False,integrationTime=30),ob.getPixelCount(row2,col2,weighted=True,integrationTime=30)
#
#        fig=plt.figure()
#        ax = fig.add_subplot(111)
#        s,bins = ob.getPixelSpectrum(row1,col1,weighted=False,integrationTime=30)
#        s2,bins = ob.getPixelSpectrum(row1,col1,weighted=True,integrationTime=30)
#
#        ax.plot(bins[0:-1],s)
#        ax.plot(bins[0:-1],s2)
#        plt.show()
#        fig=plt.figure()
#        ax = fig.add_subplot(111)
#        s,bins = ob.getPixelSpectrum(row2,col2,weighted=False,integrationTime=30)
#        s2,bins = ob.getPixelSpectrum(row2,col2,weighted=True,integrationTime=30)
#        ax.plot(bins[0:-1],s)
#        ax.plot(bins[0:-1],s2)
#        plt.show()

        startJD = ob.getFromHeader('jd')
        nSecInFile = ob.getFromHeader('exptime')
        deadMask = ob.getDeadPixels()
        for sec in range(0,nSecInFile,integrationTime):
            jd = startJD + sec/(24.*3600.)#add seconds offset to julian date
            print count,jd
            count+=1
            times.append(jd)
            titles.append('%.6f'%jd)
            #make 30s flat-fielded,hot pixel masked image
            #hotPixMaskRaw = hotPixels.findHotPixels(obsFile=ob,firstSec=sec,intTime=integrationTime)['badflag']
            #hotPixMask = hotPixels.findHotPixels(obsFile=ob,firstSec=sec,intTime=integrationTime,weighted=True)['badflag']
            frame = ob.getPixelCountImage(firstSec=sec,integrationTime=integrationTime,weighted=True)
#            hotPixMask = hotPixelsOld.findHotPixels(image=frame,firstSec=sec,intTime=integrationTime,weighted=True, display=True)['badflag']

            hotPixMask = hotPixels.checkInterval(image=frame, firstSec=sec, intTime=integrationTime, weighted=True, display=True)['mask']
#            plt.show()
#            print hotPixMask[np.isnan(hotPixMask)]
#            print hotPixMask
            
            showFrame = np.array(frame)
#            utils.plotArray(showFrame,cbar=True,normMax=np.mean(showFrame)+3*np.std(showFrame))
#            showFrame[hotPixMaskRaw == 1] = 0
#            utils.plotArray(showFrame,cbar=True)
            showFrame[hotPixMask == 1] = 0

#            utils.plotArray(showFrame,cbar=True)
            showframes.append(showFrame)

#            frame[hotPixMaskRaw == 1] = np.nan
            frame[hotPixMask == 1] = np.nan
            frame[deadMask == 0] = np.nan

#	    frame[bad_solution_mask == 1] = np.nan

            frames.append(frame)

#            print showFrame[np.isnan(showFrame)]
#            showFrame[np.isnan(showFrame)]=0
#            print showFrame[np.isnan(showFrame)]            
        
cube = np.dstack(frames)
times = np.array(times)

np.savez(npzFileName,stack=cube,jd=times)
print 'saved'
utils.makeMovie(showframes,cbar=True,frameTitles=titles,outName=gifFileName)

