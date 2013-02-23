import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
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

obsFn = FileName(run=run,date=sunsetDates[0],tstamp='20121209-120530').obs()
ob = ObsFile(obsFn)
ob.loadWvlCalFile(wvlCalFilenames[0])
ob.loadFlatCalFile(flatCalFilenames[0])
deadMask = ob.getDeadPixels()

guessX=14
guessY=8
radius=5
wvlBin=100

apertureMask=aperture(guessX,guessY,radius=radius)
bigMask = aperture(guessX,guessY,radius=radius*2)
skyMask = bigMask-apertureMask
y_values,x_values= np.where(np.logical_and(apertureMask==0,deadMask==1))
y_sky,x_sky=np.where(np.logical_and(skyMask==0,deadMask==1))

skyspectrum=[]
for i in range(len(x_sky)):
#    skyspectrum.append(ob.getPixelSpectrum(y_sky[i],x_sky[i])[0])
    skyspectrum.append(ob.getPixelSpectrum(y_sky[i],x_sky[i],wvlBinWidth=wvlBin)[0])

sky_array = np.zeros(len(skyspectrum[0]))
for j in range(len(skyspectrum[0])):
    ispectrum = np.zeros(len(skyspectrum))
    for i in range(len(skyspectrum)):    
        ispectrum[i]= skyspectrum[i][j]
    sky_array[j] = np.median(ispectrum)


#spectrum = np.zeros(len(x_values))
#bin_edges = np.zeros(len(x_values))
spectrum=[]
for i in range(len(x_values)):
#    spectrum.append(ob.getPixelSpectrum(y_values[i],x_values[i])[0]-sky_array)
    spectrum.append(ob.getPixelSpectrum(y_values[i],x_values[i],wvlBinWidth=wvlBin)[0]-sky_array)
bin_edges = ob.getPixelSpectrum(y_values[0],x_values[0],wvlBinWidth=wvlBin)[1]

summed_array = np.zeros(len(spectrum[0]))
for j in range(len(spectrum[0])):
    ispectrum = np.zeros(len(spectrum))
    for i in range(len(spectrum)):    
        ispectrum[i]= spectrum[i][j]
    summed_array[j] = np.sum(ispectrum)

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(ob.flatCalWvlBins[0:-1],summed_array)
ax.plot(bin_edges[0:-1],summed_array)
plt.xlabel('Wavelength ($\AA$)')
plt.ylabel('Counts')
plt.show()

#ob.plotPixelSpectra(8,14,weighted=True)


