#!/bin/python

'''
Author: Paul Szypryt		Date: March 8, 2013

Uses raw photon data (no integration) to calculate the Lomb-Scargle Periodogram.  
'''

import numpy as np
import scipy as sp
from scipy.signal import spectral
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
from time import time
import tables
import matplotlib.pyplot as plt
from hotpix import hotPixels


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

run = 'PAL2012'

# December 8
# First sequence, possible reflections at 12:07, 1" SE move at 12:45.
seq0 = ['120530', '121033','121536', '122039', '122542', '123045', '123548', '124051', '124554', '125057', '125601', '130103', '130606']
#seq0 = [ '121033','121536', '122039', '122542', '123045', '123548', '124051', '124554', '125057', '125601', '130103', '130606']
#seq0 = ['120530','121033']

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
seq5 = ['112709', '113212', '113714', '114216', '114718', '115220', '115722', '120224', '120727', '121229', '121732', '122234', '122736', '123238', '123740', '124242', '124744', '125246', '125748', '130250', '130752', '131254', '131756', '132258', '132800', '133303']

# Date and cal time stamp arrays
#utcDates = ['20121209','20121209','20121211', '20121211', '20121211', '20121212']
#sunsetDates = ['20121208','20121208','20121210', '20121210', '20121210', '20121211']
#calTimestamps = ['20121209-131132','20121209-133419', '20121211-074031', '20121211-074031', '20121211-133056', '20121212-133821']
utcDates = ['20121209']
sunsetDates = ['20121208']
calTimestamps = ['20121209-131132']

#seqs = [seq0,seq1,seq2,seq3,seq4,seq5]
seqs=[seq0]

timestampLists = [[utcDate+'-'+str(ts) for ts in seq] for utcDate,seq in zip(utcDates,seqs)]
wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
wvlCalFilenames[0] = '/Scratch/waveCalSolnFiles/20121210/calsol_20121211-074031.h5'
#wvlCalFilenames[1] = '/home/danica/optimusP/testing/forMatt/calsol_20121211-044853.h5'
flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'
#flatCalFilenames[1] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'

exptime = 300
app_mask = aperture(15,8,8)
y_values,x_values = np.where(app_mask==0)

timestamps =[]
tic = time()
for iSeq in range(len(seqs)):
    timestampList = timestampLists[iSeq]
    wfn = wvlCalFilenames[iSeq]
    ffn = flatCalFilenames[iSeq]
    sunsetDate = sunsetDates[iSeq]
    for i,ts in enumerate(timestampList):
        print 'Loading',ts
        obsFn = FileName(run=run,date=sunsetDate,tstamp=ts).obs()
        ob = ObsFile(obsFn)
        ob.loadWvlCalFile(wfn)
        ob.loadFlatCalFile(ffn)
	ob.setWvlCutoffs(3000,5000)
	for j in range(len(x_values)):
	    x=ob.getPixelWvlList(iRow=y_values[j],iCol=x_values[j])
	    timestamps = np.append(timestamps,i*exptime + (x['timestamps']))
print 'Loaded files in',time()-tic, 'seconds.'

# Put photon in events in small time bins ~ 1 microsecond goal
tic = time()
print 'Binning photon events...'
total_seconds = len(seq0*exptime)
bins_per_second = 10**2
#time_bins = np.linspace(0,total_seconds,total_seconds*bins_per_second+1)
#jd = time_bins/86400
#counts_per_timestep = np.zeros(len(time_bins)-1)
#for j in range(len(time_bins)-1):
    #counts_per_timestep[j]=len(np.where(np.logical_and(timestamps>=time_bins[j],timestamps<time_bins[j+1]))[:][0])
counts_per_timestep, bin_edges = np.histogram(timestamps,bins = bins_per_second*total_seconds)
jd =bin_edges[0:-1]/86400
print 'Finished binning in',time()-tic, 'seconds.'
 
scaled_counts = (counts_per_timestep-counts_per_timestep.mean())/counts_per_timestep.std()

# Create array of frequencies to check for Fourier components, function requires angular frequencies
freqs=np.linspace(0.1,250,10000)
angular_freqs=2*np.pi*freqs
print 'Calculating Fourier components...'
tic = time()
periodogram = spectral.lombscargle(jd, scaled_counts, angular_freqs)
print 'Calculated Fourier components in',time()-tic, 'seconds.'

# Calculate eclipse period and frequency, and compare it to expected value
eclipse_period = 2*np.pi/(angular_freqs[np.argmax(periodogram)])
eclipse_frequency = 1/eclipse_period
expected_period = 0.01966127 # in days
print 'Eclipse period =',eclipse_period,'days.'
print 'Eclipse frequency =',eclipse_frequency, 'cycles/day.'
print 'Percent error = ' + str(100*(eclipse_period-expected_period)/expected_period) + '%'

# Create a figure with light curve in top plot and periodogram in bottom plot
fig = plt.figure()
# Plot light curve
ax = fig.add_subplot(211)
#ax.plot(jd,scaled_counts,'b.')
ax.set_title('Light Curve')
plt.xlabel('Time (Days)')
plt.ylabel('Scaled Counts')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
# Plot fourier transform
ax = fig.add_subplot(212)
ax.plot(freqs,periodogram)
ax.set_title('Periodogram')
plt.xlabel('Frequency (Cycles/Day)')
plt.ylabel('Transform Component')
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] + ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(16)
plt.show()





