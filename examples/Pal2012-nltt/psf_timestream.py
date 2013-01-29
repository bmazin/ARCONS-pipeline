#!/bin/python

import numpy as np
from util.ObsFile import ObsFile
from util.FileName import FileName
from util import utils
import tables
import matplotlib.pyplot as plt


class ApertureList(tables.IsDescription):
    ArrivalTime = tables.Float64Col()
    Flag = UInt8Col()
    Phase = tables.Float32Col()
    PixelRow = tables.UInt8Col()
    PixelCol = tables.UInt8Col()
    Wavelength = tables.Float32Col()
    Weight = tables.Float32Col()

def circ(startpx,startpy,radius=3):
    r = radius
    length = 2*r 
    height = length
    allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
    ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
    pixx = []
    pixy = []
    for x in allx:
        for y in ally:
            if (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 <= (r)**2:
                pixx.append(x)
                pixy.append(y)
    return pixx,pixy

run = 'PAL2012'
firstRun20121210 = ['022500','023005','023507','024009','024511','025013','025515','030017','030519','031021']
secondRun20121210 = ['031755','032257','032759','033301','033803','034305','034807','035309','035811','040313','040815','041318','041820','042322','042824','043326','043828','044331']
firstRun20121211 = ['023841','024343','024845','025348','025850','030352','030854','031356','031858']

utcDate = '20121211'
sunsetDate = '20121210'
#timestampListPost = firstRun20121210
#timestampListPost = secondRun20121210
timestampListPost=firstRun20121211
#circCol,circRow = circ(28,32,6)#firstRun20121210
#circCol,circRow = circ(29,31,6)#secondRun20121210
circCol,circRow = circ(14,7,6)#firstRun20121211
runlabel = '1'
out = '/home/mstrader/data/nltt/'+'timestream'+runlabel+'.txt'
outh5 = '/home/mstrader/data/nltt/'+'list'+runlabel+'.h5'


timestampList = [utcDate+'-'+ts for ts in timestampListPost]

files = []
pl = []
times = []
for i,ts in enumerate(timestampList):
    print 'loading',ts
    fn = FileName(run=run,date=sunsetDate,tstamp=ts).photonList()
    obsFn = FileName(run=run,date=sunsetDate,tstamp=ts).obs()
    files.append(tables.openFile(fn,mode='r'))
    ob = ObsFile(obsFn)
    ut = ob.getFromHeader('unixtime')
    pl.append(files[i].root.photons.photons.read())
    times.append(np.array(pl[i]['ArrivalTime'],dtype=np.float64) + ut)


all=np.concatenate(pl)
allTimes = np.concatenate(times)

nlttPSFByPixels = []
nlttTimesByPixels = []
img = np.zeros((46,44))
rawImg = np.zeros((46,44))
for iPixel in range(len(circCol)):
    x = circCol[iPixel]
    y = circRow[iPixel]
    timesInPixel = allTimes[np.logical_and(all['Xpix']==x,all['Ypix']==y)]
    inPixel = all[np.logical_and(all['Xpix']==x,all['Ypix']==y)]
    nlttPSFByPixels.append(inPixel)
    nlttTimesByPixels.append(timesInPixel)
    img[y,x] = sum(inPixel['FlatWeight'])
    rawImg[y,x] = len(inPixel)

nlttPSF = np.concatenate(nlttPSFByPixels)
nlttPSFTimes = np.concatenate(nlttTimesByPixels)
secsInDay = 24*60*60.
period = 0.2350606*secsInDay
phases = (nlttPSFTimes % period)/(period)
newtype=[('ArrivalTime', '<f8'),('Flag', '|u1'), ('Phase', '<f4'), ('PixelCol', '|u1'), ('PixelRow', '|u1'), ('Wavelength', '<f4'), ('Weight', '<f4')]
nPhotons = len(nlttPSF)
newTable = np.recarray(nPhotons,dtype=newtype)
newTable['ArrivalTime'] = nlttPSFTimes
newTable['Flag'] = nlttPSF['Flag']
newTable['Phase'] = phases



timestream,timeEdges = np.histogram(nlttPSFTimes,weights=nlttPSF['FlatWeight'],bins=300*len(timestampList))
phaseTimestream,phaseTimeEdges = np.histogram(phases,weights=nlttPSF['FlatWeight'],bins=30*len(timestampList))

plt.plot(timeEdges[:-1],timestream)
plt.show()
counts = [len(pixelPhotons) for pixelPhotons in nlttPSFByPixels]
utils.plotArray(img,cbar=True,normMax=800000)
utils.plotArray(rawImg,cbar=True)
print img[31,29]
print rawImg[31,29]
tbl = np.vstack([timeEdges[:-1],timestream])
np.savetxt(out,tbl.T)
