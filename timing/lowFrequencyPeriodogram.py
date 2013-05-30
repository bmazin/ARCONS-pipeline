#!/bin/python

'''
Author: Paul Szypryt		Date: March 8, 2013

Uses raw photon data (no fits, too small timescales) to calculate the Lomb-Scargle Periodogram.  Does this in a specified
aperture of given x and y location and radius.  Applies wavelength cals, flat cals, and other typical ObsFiles functions.
'''

import numpy as np
import scipy as sp
import sys
from scipy.signal import spectral
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
from time import time
import tables
import os
import matplotlib.pyplot as plt
import hotpix.hotPixels as hp


# Function that creates an aperture mask of specified radius at position startpx, startpy.
def aperture(startpx,startpy,radius=3):
    r = radius
    length = 2*r 
    height = length
    allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
    ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
    pixx = []
    pixy = []
    mask=np.ones((46,44))
    for x in allx:
        for y in ally:
            if (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 <= (r)**2 and 0 <= y and y < 46 and 0 <= x and x < 44:
                mask[y,x]=0.
    return mask

# Name of the observing run.
run = 'PAL2012'

# Name of save path.
savePath= os.getenv("HOME") + '/Scratch/'

# SDSS J0926 data, data sequences for other objects can be similarly created.
# December 8
# First sequence, possible reflections at 12:07, 1" SE move at 12:45.
seq0 = ['120530', '121033','121536', '122039', '122542', '123045', '123548', '124051', '124554', '125057', '125601', '130103', '130606']
# Sequence during warming up, may need to omit.
seq1 = ['131254', '131802', '132304', '132807']
# December 10
# Sequence during 2nd day on object. Moved to confirm position in 082422 between 90 and 105s.'080916' hot pix.
seq2 = ['074405', '074907', '075410', '075912', '080414', '080916', '081418', '081920', '082422']
# Refocused and started guiding again at 8:37. Ommitting seq 083451, which occured during refocus.
seq3 = ['084029', '084532', '085034', '085536', '090038']
# Back toward end of night, not sure about whether to use cal file '20121211-115429' or '20121211-133056'.  Also, may need to cut out last 2 obs files.
seq4 = ['120152', '120654', '121157', '121700', '122203', '122706', '123209', '123712', '124215', '124809', '125312', '125814', '130316', '130818', '131320', '131822', '132324']
# December 11
# Final sequence, toward end of run, thin high clouds at around 12:50, moved to confirm position at '122234', also at '130752' at 125s.
seq5 = ['112709', '113212', '113714', '114216', '114718', '115220', '115722', '120224', '120727', '121229', '121732', '122234', '122736', '123238', '123740', '124242', '124744', '125246', '125748', '130250', '130752', '131254', '131756', '132258', '132800', '133303']

# Specify sunset dates (date of beginning of an observation, which also corresponds to folder names of data products) and utc dates.
sunsetDates = ['20121208']
utcDates = ['20121209']

# Specify which wavelength calibration file to use.
calTimestamps = ['20121209-131132']

# Pick out which of the defined sequences to use in the periodogram.
seqs=[seq0]

# Append observation time in sequence to utc date.  Used to create the observation file names.
timestampLists = [[utcDate+'-'+str(ts) for ts in seq] for utcDate,seq in zip(utcDates,seqs)]

# Create wavelength and flat cal file names
wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
# No twilights taken on December 8, using December 7 twilights to do flat cal instead.
flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

# Specify the amount of time bins to make per second of observations.  E.g. binsPerSecond = 10**6 corresponds to microsecond binning.
binsPerSecond = 0.1

# Specify frequencies for which the periodogram algorithm should find transform components.
freqs = np.linspace(10**-4,10**-3,num=10**4)
angularFreqs=2*np.pi*freqs

# Define the area of an observation containing the object.
appMask = aperture(15,8,8)
yValues,xValues = np.where(appMask==0)

# Initialize some variables.
unixOffset=0.0
countsPerTimestep=[]
timeArray=[]
tic = time()

for iSeq in range(len(seqs)):
    # Pick out timestamps, wave cal, flat cals, and sunset date for a particular sequence.
    timestampList = timestampLists[iSeq]
    wfn = wvlCalFilenames[iSeq]
    ffn = flatCalFilenames[iSeq]
    sunsetDate = sunsetDates[iSeq]
    for i,ts in enumerate(timestampList):
        print 'Loading',ts
        timestamps =[]
        # Create ObsFile instance.
        obsFn = FileName(run=run,date=sunsetDate,tstamp=ts).obs()
        ob = ObsFile(obsFn)
        # Load roach time delay corrections.
        ob.loadTimeAdjustmentFile(FileName(run=run).timeAdjustments())
        # Retrieve exposure time and unix time from obs file header.
        exptime = ob.getFromHeader('exptime')
        unixtime= ob.getFromHeader('unixtime')
        # Set unix time of first obs file as the zero point.
        if i == 0 and iSeq == 0:
            unixOffset = unixtime
        # Search for time mask for given observation file.  If the time mask does not exist, create it.
        index1 = obsFn.find('_')
        hotPixFn = '/Scratch/timeMasks/timeMask' + obsFn[index1:]
        if not os.path.exists(hotPixFn):
            hp.findHotPixels(obsFn,hotPixFn)
            print "Flux file pixel mask saved to %s"%(hotPixFn)
        # Load time mask, wavelength calibration, and flat calibration and set wavelenth cutoffs.
        ob.loadHotPixCalFile(hotPixFn,switchOnMask=True)
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
        ob.setWvlCutoffs(3000,8000)
        # Get a list of timestamps with photons from each pixel in the aperture.
        for j in range(len(xValues)):
            x=ob.getPixelWvlList(iRow=yValues[j],iCol=xValues[j])
            timestamps = np.append(timestamps,x['timestamps'])
        # Histogram the timestamps with bins corresponding to the binsPerSecond variable.
        binnedCounts, binEdges = np.histogram(timestamps,range=[0,exptime],bins = binsPerSecond*exptime)
        # Append photon count data to a total list across multiple obs files.
        countsPerTimestep = np.append(countsPerTimestep,binnedCounts)
        # Use the unix time to determine the delay between observations and create an accurately spaced time axis.
        binEdges=binEdges[0:-1]
        timeArray=np.append(timeArray, (unixtime-unixOffset) + binEdges)
	   			
print 'Loaded files in',time()-tic, 'seconds.'

# Normalize the photon count data.  This is required for the periodogram algorithm.
scaledCounts = (countsPerTimestep-countsPerTimestep.mean())/countsPerTimestep.std()

# Use timeArray, scaledCounts, and angularFreqs to create a Lomb-Scargle periodogram.
print 'Calculating Fourier components...'
tic = time()
periodogram = spectral.lombscargle(timeArray, scaledCounts, angularFreqs)
print 'Calculated Fourier components in',time()-tic, 'seconds.'

# Save data to txt files.
f=open(savePath+'timeDataTest.txt','w')
g=open(savePath+'frequencyDataTest.txt','w')
for time in range(len(timeArray)):
    f=open(savePath + 'timeDataTest.txt','a')
    f.write(str(timeArray[time]) + '\t' + str(scaledCounts[time]) + '\n')
for iFreq in range(len(freqs)):
    g=open(savePath + 'frequencyDataTest.txt','a')
    g.write(str(freqs[iFreq]) + '\t' + str(periodogram[iFreq]) + '\n')
f.close()
g.close()

fig = plt.figure()
# Plot light curve
ax = fig.add_subplot(211)
ax.plot(timeArray,scaledCounts,'b.')
ax.set_title('Light Curve')
plt.xlabel('Julian Date')
plt.ylabel('Scaled Counts')
# Plot fourier transform
ax = fig.add_subplot(212)
ax.plot(freqs,periodogram)
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_title('Periodogram')
plt.xlabel('Frequency (Cycles/Day)')
plt.ylabel('Power')
plt.show()



