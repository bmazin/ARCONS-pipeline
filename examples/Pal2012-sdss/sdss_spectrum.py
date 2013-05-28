import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
import os
import hotpix.hotPixels as hp
import tables
import matplotlib.pyplot as plt

def getTimeMaskFileName(obsFileName):
    scratchDir = os.getenv('INTERM_PATH')
    hotPixDir = os.path.join(scratchDir,'timeMasks')
    fileName = obsFileName.split('/')[-1]
    fileNameBase = fileName.split('_')[-1]
    newName = 'timeMask_'+fileNameBase
    fullHotPixFileName = os.path.join(hotPixDir,newName)
    return fullHotPixFileName

run = 'PAL2012'

# December 8
# First sequence, possible reflections at 12:07, 1" SE move at 12:45.
seq0 = ['120530', '121033','121536', '122039', '122542', '123045', '123548', '124051', '124554', '125057', '125601', '130103', '130606']

# Sequence during warming up, may need to omit.
seq1 = ['131254', '131802', '132304', '132807']

seqs = [seq0,seq1]

#calTimestamps = ['20121209-131132','20121209-133419', '20121211-074031', '20121211-074031', '20121211-133056', '20121212-133821']
#utcDates = ['20121209','20121209','20121211', '20121211', '20121211', '20121212']
#sunsetDates = ['20121208','20121208','20121210', '20121210', '20121210', '20121211']

calTimestamps = ['20121209-131132','20121209-133419']
utcDates = ['20121209','20121209']
sunsetDates = ['20121208','20121208']

timestampLists = [[utcDate+'-'+str(ts) for ts in seq] for utcDate,seq in zip(utcDates,seqs)]
wvlCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).calSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]
flatCalFilenames = [FileName(run=run,date=sunsetDate,tstamp=calTimestamp).flatSoln() for sunsetDate,calTimestamp in zip(sunsetDates,calTimestamps)]

#wvlCalFilenames[0] = '/Scratch/waveCalSolnFiles/20121210/calsol_20121211-074031.h5'
#wvlCalFilenames[1] = '/home/danica/optimusP/testing/forMatt/calsol_20121211-044853.h5'
flatCalFilenames[0] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'
flatCalFilenames[1] = '/Scratch/flatCalSolnFiles/20121207/flatsol_20121207.h5'
fluxCalFileNames = ['/Scratch/fluxCalSolnFiles/20121206/fluxsol_20121207-124034.h5']

obsFn = FileName(run=run,date=sunsetDates[0],tstamp='20121209-120530').obs()
ob = ObsFile(obsFn)
print 'Loading wavelength calibration solution: ' + wvlCalFilenames[0]
ob.loadWvlCalFile(wvlCalFilenames[0])
print 'Loading flat calibration solution: ' + flatCalFilenames[0]
ob.loadFlatCalFile(flatCalFilenames[0])
ob.loadFluxCalFile(fluxCalFileNames[0])

#load/generate hot pixel mask file
HotPixFile = getTimeMaskFileName(obsFn)
if not os.path.exists(HotPixFile):
    hp.findHotPixels(obsFn,HotPixFile)
    print "Flux file pixel mask saved to %s"%(HotPixFile)
ob.loadHotPixCalFile(HotPixFile)
print "Hot pixel mask loaded %s"%(HotPixFile)

frame = ob.getPixelCountImage(firstSec=0,integrationTime=300,weighted=True)
#hotPixMask = hotPixels.checkInterval(image=frame, firstSec=0, intTime=300, weighted=True, display=False)['mask']

#summed_array,bin_edges=ob.getApertureSpectrum(pixelCol=14,pixelRow=8,radius=7)
ob.plotApertureSpectrum(pixelCol=14,pixelRow=8,radius=7,weighted = True,fluxWeighted=True,lowCut=3000,highCut=9000)


'''
h = 6.626068*10**-34
c = 299792458.0
k = 1.3806503*10**-23
T = 10000.0

numerator = 6*10**-29
denominator= (bin_edges*10**-10)**5*(np.exp(((h*c)/(bin_edges*k*T*10**-10))) - 1)
bbcurve = numerator/denominator
#print bbcurve


fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(ob.flatCalWvlBins[0:-1],summed_array)
#ax.plot(bin_edges[0:-1],summed_array,bin_edges[0:-1],bbcurve[0:-1])
ax.plot(bin_edges[12:-2],summed_array[12:-1])
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Counts')
plt.show()
'''

#ob.plotPixelSpectra(8,14,weighted=True)


