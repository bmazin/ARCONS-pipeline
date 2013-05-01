#!/bin/python

'''
Author: Paul Szypryt		Date: April 19, 2013

Uses raw photon data (no fits, too small timescales) to calculate the Lomb-Scargle Periodogram.  Does this in a specified
aperture of given x and y location and radius.  Applies wavelength cals, flat cals, and other typical ObsFiles functions.
Averages over selected timesteps, optimized for high frequencies (>1 Hz).
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

'''
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
'''

# Select which obs and calibration files to use.
# Specify obs file timestamp to use in periodogram.
obsTimestamp = '20121209-121033'
# Specify sunset dates (date of beginning of an observation, which also corresponds to folder names of data products) and utc dates.
sunsetDate = '20121208'
utcDate = '20121209'
# Specify which wavelength calibration file to use.
calTimestamp = '20121209-131132'
# Create wavelength and flat cal file names
wvlCalFilename = FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln()
flatCalFilename = FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln()
# No twilights taken on December 8, using December 7 twilights to do flat cal instead.
flatCalFilename = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

# Run standard ObsFile functions.
# Create ObsFile instance.
tic = time()
print 'Loading obs file and performing calibrations ...'
obsFn = FileName(run=run,date=sunsetDate,tstamp=obsTimestamp).obs()
ob = ObsFile(obsFn)
# Load roach time delay corrections.
ob.loadTimeAdjustmentFile(FileName(run=run).timeAdjustments())
# Search for time mask for given observation file.  If the time mask does not exist, create it.
index1 = obsFn.find('_')
hotPixFn = '/Scratch/timeMasks/timeMask' + obsFn[index1:]
if not os.path.exists(hotPixFn):
    hp.findHotPixels(obsFn,hotPixFn)
    print "Flux file pixel mask saved to %s"%(hotPixFn)        
# Load time mask, wavelength calibration, and flat calibration and set wavelenth cutoffs.
ob.loadHotPixCalFile(hotPixFn,switchOnMask=True)
ob.loadWvlCalFile(wvlCalFilename)
ob.loadFlatCalFile(flatCalFilename)
ob.setWvlCutoffs(3000,5000)
print 'Total load time: ' + str(time()-tic) + 's'

tic = time()
print 'Retrieving photon timestamps...'
# Define the area of an observation containing the object.
appMask = aperture(15,8,8)
yValues,xValues = np.where(appMask==0)
# Get a list of timestamps with photons from each pixel in the aperture.
timestamps=[]
for j in range(len(xValues)):
    x=ob.getPixelWvlList(iRow=yValues[j],iCol=xValues[j])
    timestamps = np.append(timestamps,x['timestamps'])
print 'Total retrieval time: ' + str(time()-tic) + 's'

# Use timestamp data to perform periodograms.
timestep = 5.0*10**-4
averagingTime = 2.0
exptime = ob.getFromHeader('exptime')
totalAverageBins = int(exptime/averagingTime)
timestepsPerBin = int(averagingTime/timestep)
# Specify frequencies for which the periodogram algorithm should find transform components.
freqs = np.linspace(1,1000,num=10**4)
angularFreqs=2*np.pi*freqs

binnedCounts, binEdges = np.histogram(timestamps,range=[0,exptime],bins = exptime/timestep)
times = binEdges[0:timestepsPerBin]

periodogramStack=np.zeros((totalAverageBins,len(freqs)))
tic = time()
print 'Calculating individual Fourier Transforms...'
for bin in range(int(totalAverageBins)):
    counts = binnedCounts[bin*timestepsPerBin:(bin+1)*timestepsPerBin]
    scaledCounts = (counts-counts.mean())/counts.std()
    periodogram = spectral.lombscargle(times, scaledCounts, angularFreqs)
    periodogramStack[bin,:] = periodogram
    completed = 100*(float(bin)+1)/float(totalAverageBins)
    print '%.1f' % completed + '% complete'
print 'Total Fourier Transform time: ' + str(time()-tic) + 's'

averageStack = np.zeros(len(freqs))
for i in range(len(freqs)):
    timeInterval=np.zeros(totalAverageBins)
    for j in range(totalAverageBins):
	    timeInterval[j]=periodogramStack[j][i]
    averageStack[i]=np.average(timeInterval)
	   			
print 'Saving data to txt files...'
# Save data to txt files.
g=open(savePath+'frequencyData.txt','w')
for iFreq in range(len(freqs)):
    g=open(savePath + 'frequencyData.txt','a')
    g.write(str(freqs[iFreq]) + '\t' + str(averageStack[iFreq]) + '\n')
g.close()

# Plot the periodogram transform components (power?) vs frequencies.
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(freqs,averageStack)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_title('Periodogram')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Transform Component')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
plt.show()



