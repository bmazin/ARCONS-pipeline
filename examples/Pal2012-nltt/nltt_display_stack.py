#!/bin/python

import numpy as np
from util.ObsFile import ObsFile
from util.FileName import FileName
from util import utils
import tables
import matplotlib.pyplot as plt
from hotpix import hotPixels



run = 'PAL2012'
seq0 = ['023005','023507','024009','024511','025013','025515','030017','030519','031021']
seq1 = ['031755','032257','032759','033301','033803','034305','034807','035309','035811','040313','040815','041318','041820','042322','042824','043326','043828','044331']
seq2 = ['023841','024343','024845','025348','025850','030352','030854','031356','031858']

seqs = [seq0,seq1,seq2]
utcDates = ['20121211','20121211','20121212']
sunsetDates = ['20121210','20121210','20121211']
timestampLists = [[utcDate+'-'+ts for ts in seq] for utcDate,seq in zip(utcDates,seqs)]

calTimestamps = ['20121211-020441','20121211-044853','20121212-032455']
wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
wvlCalFilenames[0] = '/home/danica/optimusP/testing/forMatt/calsol_20121211-020441.h5'
wvlCalFilenames[1] = '/home/danica/optimusP/testing/forMatt/calsol_20121211-044853.h5'
flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]

integrationTime=30
frames = []
showframes = []
times = []
titles = []
#plt.ion()
count= 0
print (len(seq0)+len(seq1)+len(seq2))*300/integrationTime,'frames to make'
for iSeq in range(3):
    timestampList = timestampLists[iSeq]
    wfn = wvlCalFilenames[iSeq]
    ffn = flatCalFilenames[iSeq]
    sunsetDate = sunsetDates[iSeq]
    for i,ts in enumerate(timestampList):
        print 'loading',ts
        obsFn = FileName(run=run,date=sunsetDate,tstamp=ts).obs()
        ob = ObsFile(obsFn)
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
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
            hotPixMask = hotPixels.findHotPixels(image=frame,firstSec=sec,intTime=integrationTime,weighted=True)['badflag']
#            plt.show()

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
            frames.append(frame)
            
        
cube = np.dstack(frames)
times = np.array(times)
np.savez('nlttImageStack6000.npz',stack=cube,jd=times)
utils.makeMovie(showframes,cbar=True,frameTitles=titles,outName='nlttImageStack6000.gif')

