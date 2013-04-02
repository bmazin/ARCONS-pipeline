#!/bin/python

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
import tables
import matplotlib.pyplot as plt
import hotpix.hotPixels as hp
import os
from time import time


run = 'PAL2012'

# December 8
# First sequence, possible reflections at 12:07, 1" SE move at 12:45. (14,8)
seq0 = ['120530', '121033','121536', '122039', '122542', '123045', '123548', '124051', '124554', '125057', '125601', '130103', '130606']

# Sequence during warming up, may need to omit. (14,8)
seq1 = ['131254', '131802', '132304', '132807']

# December 10
# Sequence during 2nd day on object. Moved to confirm position in 082422 between 90 and 105s.'080916' hot pix (31,30)
seq2 = ['074405', '074907', '075410', '075912', '080414', '080916', '081418', '081920', '082422']

# Refocused and started guiding again at 8:37. Ommitting seq 083451, which occured during refocus. (31,30)
seq3 = ['084029', '084532', '085034', '085536', '090038']

# Back toward end of night, not sure about whether to use cal file '20121211-115429' or '20121211-133056'.  Also, may need to cut out last 2 obs files (14,8)
seq4 = ['120152', '120654', '121157', '121700', '122203', '122706', '123209', '123712', '124215', '124809', '125312', '125814', '130316', '130818', '131320', '131822', '132324']

# December 11
# Final sequence, toward end of run, thin high clouds at around 12:50, moved to confirm position at '122234', also at '130752' at 125s. (16,15)
seq5 = ['112709', '113212', '113714', '114216', '114718', '115220', '115722', '120224', '120727', '121229', '121732', '122234', '122736', '123238', '123740', '124242', '124744', '125246', '125748', '130250', '130752', '131254', '131756', '132258', '132800', '133303']

seqHMCnc= ['074719', '075225','075753','080456','080959','081501', '082004', '082507', '083009','083747', '084257', '084759', '085302', '085805','090308','090810']
seq0651=['095930','100453','100955','101457','101959', '102501', '103003', '103505', '104007','104509', '105050', '105553', '110055','110557','111100','111602', '112104','112606', '113108', '113611', '114113', '114615', '115118']

# Date and cal time stamp arrays
#utcDates = ['20121209','20121209','20121211', '20121211', '20121211', '20121212']
#sunsetDates = ['20121208','20121208','20121210', '20121210', '20121210', '20121211']
#calTimestamps = ['20121209-131132','20121209-133419', '20121211-074031', '20121211-074031', '20121211-133056', '20121212-133821']
utcDates = ['20121209','20121209']
sunsetDates = ['20121208','20121208']
calTimestamps = ['20121209-131132','20121209-133419']
#utcDates = ['20121206']
#sunsetDates = ['20121205']
#calTimestamps = ['20121206-115709']

#seqs = [seq0,seq1,seq2,seq3,seq4,seq5]
seqs=[seq0,seq1]

timestampLists = [[utcDate+'-'+str(ts) for ts in seq] for utcDate,seq in zip(utcDates,seqs)]
wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
#wvlCalFilenames[0] = '/Scratch/waveCalSolnFiles/20121210/calsol_20121211-074031.h5'
#wvlCalFilenames[1] = '/home/danica/optimusP/testing/forMatt/calsol_20121211-044853.h5'
flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
#flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121206/flatsol_20121206.h5'
flatCalFilenames[1] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

#/Scratch/waveCalSolnFiles/20121208/calsol_20121209-131132.h5

integrationTime=10
frames = []
showframes = []
times = []
titles = []
expTime = 300
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
	index1 = obsFn.find('_')
	hotPixFn = '/Scratch/timeMasks/timeMask' + obsFn[index1:]
        if not os.path.exists(hotPixFn):
            hp.findHotPixels(obsFn,hotPixFn)
            print "Flux file pixel mask saved to %s"%(hotPixFn)
        ob.loadHotPixCalFile(hotPixFn,switchOnMask=True)
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
        ob.setWvlCutoffs(3000,5000)

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
        deadMask = ob.getDeadPixels()
	#print 'Dead mask load time = ', time()-tic
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
np.savez('/home/pszypryt/sdss_data/test.npz',stack=cube,jd=times)
utils.makeMovie(showframes,cbar=True,frameTitles=titles,outName='/home/pszypryt/sdss_data/test.gif')

