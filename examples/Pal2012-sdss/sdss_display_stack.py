#!/bin/python

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
import tables
import matplotlib.pyplot as plt
from hotpix import hotPixels


run = 'PAL2012'

# December 8
# First sequence, possible reflections at 12:07, 1" SE move at 12:45.
seq0 = ['120530', '121033','121536', '122039', '122542', '123045', '123548', '124051', '124554', '125057', '125601', '130103', '130606']

# Sequence during warming up, may need to omit.
seq1 = ['131254', '131802', '132304', '132807']

# December 10
# Sequence during 2nd day on object. Moved to confirm position in 082422 between 90 and 105s.'080916' hot pix
seq2 = ['074405', '074907', '075410', '075912', '080414', '080916', '081418', '081920', '082422']

# Refocused and started guiding again at 8:37. Ommitting seq 083451, which occured during refocus.
seq3 = ['084029', '084532', '085034', '085536', '090038']

# Back toward end of night, not sure about whether to use cal file '20121211-115429' or '20121211-133056'.  Also, may need to cut out last 2 obs files
seq4 = ['120152', '120654', '121157', '121700', '122203', '122706', '123209', '123712', '124215', '124809', '125312', '125814', '130316', '130818', '131320', '131822', '132324']

# December 11
# Final sequence, toward end of run, thin high clouds at around 12:50, moved to confirm position at '122234', also at '130752' at 125s.
seq5 = ['112709', '113212', '113714', '114216', '114718', '115220', '115722', '120224', '120727', '121229', '121732', '122234', '122736', '123238']
seq6 = ['123740', '124242', '124744', '125246', '125748', '130250', '130752', '131254', '131756', '132258', '132800', '133303']

# Date and cal time stamp arrays
#utcDates = ['20121209', '20121209', '20121211', '20121211', '20121211', '20121212']
#sunsetDates = ['20121208', '20121208', '20121210', '20121210', '20121210', '20121211']

#calTimestamps = ['20121209-131132','20121209-133419', '20121211-090613', '20121211-090613', '20121211-133056', '20121212-111847'] 
#'20121212-133821' is terrible! 
#try replacing '20121211-074031' with '20121211-090613'
utcDates = ['20121211', '20121211']
sunsetDates = ['20121210', '20121210']
#utcDates = ['20121209', '20121209', '20121211', '20121211', '20121211', '20121212']
#sunsetDates = ['20121208', '20121208', '20121210', '20121210', '20121210', '20121211']

calTimestamps = ['20121211-090613', '20121211-090613']
seqs = [seq2,seq3]

timestampLists = [[utcDate+'-'+str(ts) for ts in seq] for utcDate,seq in zip(utcDates,seqs)]

wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
#wvlCalFilenames[0] = '/Scratch/waveCalSolnFiles/20121210/calsol_20121211-074031.h5'
#wvlCalFilenames[1] = '/home/danica/optimusP/testing/forMatt/calsol_20121211-044853.h5'
#flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(['20121210','20121210'],['20121211-074031','20121211-074031'])]
flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
#flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'
#flatCalFilenames[1] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

#/Scratch/waveCalSolnFiles/20121208/calsol_20121209-131132.h5

integrationTime=3
frames = []
showframes = []
times = []
titles = []
expTime = 300
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
    wfn = wvlCalFilenames[iSeq]
    ffn = flatCalFilenames[iSeq]
    sunsetDate = sunsetDates[iSeq]
    for i,ts in enumerate(timestampList):
        print 'loading',ts
        obsFn = FileName(run=run,date=sunsetDate,tstamp=ts).obs()
        ob = ObsFile(obsFn)
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
        ob.setWvlCutoffs(5000,7000)

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
            frames.append(frame)

#            print showFrame[np.isnan(showFrame)]
#            showFrame[np.isnan(showFrame)]=0
#            print showFrame[np.isnan(showFrame)]            
        
cube = np.dstack(frames)
times = np.array(times)
np.savez('/Scratch/dataProcessing/SDSS_J0926/AllData/FirstDec10SIImageStackRednewCal.npz',stack=cube,jd=times)
print 'saved'
utils.makeMovie(showframes,cbar=True,frameTitles=titles,outName='/Scratch/dataProcessing/SDSS_J0926/AllData/FirstDec10SIImageStackRednewCal.gif')

