import numpy as np
from util.FileName import FileName
from util import utils
from util.popup import plotArray,pop
import scipy.stats
import tables
import os
import re

#  8 bits - channel
# 12 bits - Parabola Fit Peak Height
# 12 bits - Sampled Peak Height
# 12 bits - Low pass filter baseline
# 20 bits - Microsecond timestamp

nBitsAfterParabolaPeak = 44
nBitsAfterBaseline = 20
nBitsInPulseHeight = 12
nBitsInTimestamp = 20

#bitmask of 12 ones
pulseMask = int(nBitsInPulseHeight * '1', 2)
#bitmask of 20 ones
timestampMask = int(nBitsInTimestamp * '1', 2)
tickDuration = 1.e-6 #s
deleteBaseMask = int((8+12+12)*'1'+12*'0'+20*'1',2)

def parsePhotonPackets(packets):
        """
        Parses an array of uint64 packets with the obs file format
        inter is an interval of time values to mask out
        returns a list of timestamps,parabolaFitPeaks,baselines
        """
        # parse all packets
        packetsAll = np.array(packets, dtype='uint64') #64 bit photon packet
        timestamps = np.bitwise_and(packets, timestampMask)
        baselines = np.bitwise_and(np.right_shift(packets, nBitsAfterBaseline), pulseMask)
        return {'timestamps':timestamps, 'baselines':baselines}

def getTimedPacketList(pixelData):
    timestamps = []
    baselines = []
    nRows = len(pixelData)
    rowBreaks = np.zeros(nRows-1)
    breakIndex = 0
    for t in range(nRows):
            parseDict = parsePhotonPackets(pixelData[t])
            times = parseDict['timestamps']
            bases = parseDict['baselines']
            times = tickDuration * times + t

            #find where the final array will have to be split to make a list of arrays with the shape of pixelData
            breakIndex += len(times)
            if t < nRows-1:
                rowBreaks[t] = breakIndex
            timestamps.append(times)
            baselines.append(bases)
    if len(pixelData) > 0:
        timestamps = np.concatenate(timestamps)
        baselines = np.concatenate(baselines)
    else:
        timestamps = np.array([])
        baselines = np.array([])

    return {'timestamps':timestamps,'baselines':baselines,'rowBreakIndices':rowBreaks}

def smoothBaselines(baselines,timestamps,timeToAverage=500e-3,upperTrimFraction=0.5):
    #A time based smoothing
    #For each baseline point, all points within a given time interval around the point is
    #trimmed and averaged to make a new baseline point
    modBases = np.array(baselines,dtype=np.uint32)
    for i,ts in enumerate(timestamps):
        startIdx = np.searchsorted(timestamps,ts-timeToAverage/2.)
        endIdx = np.searchsorted(timestamps,ts+timeToAverage/2.)
        trimmedSample = scipy.stats.mstats.trim(baselines[startIdx:endIdx],(upperTrimFraction,0),relative=True)
        if np.all(trimmedSample.mask) == False: #if there are some points that are not trimmed, take the average
            modBases[i] = int(np.ma.mean(trimmedSample))
    return modBases

def smoothObs(run,sunsetDate,timestamp):
    smoothRun = run+'s'
    obsFileName = FileName(run=smoothRun,date=sunsetDate,tstamp=timestamp).obs()
    obs = tables.openFile(obsFileName,mode='a')
    #get the name of the datasets
    for group in obs.root:
        name = group._v_name
        matchObj = re.match(r'r\d+',name)
        if matchObj != None: #This is a roach group
            roachGroup = group
            roachNum = int(name[1:])
            print 'roach',roachNum
            for pixelGroup in roachGroup:
                print pixelGroup._v_name,
                for dataset in pixelGroup: #There is only one
                    print len(dataset),
                    data = dataset.read()
                    parsedDict = getTimedPacketList(data)
                    timestamps = parsedDict['timestamps']
                    baselines = parsedDict['baselines']
                    print len(timestamps)
                    modBases = smoothBaselines(baselines=baselines,timestamps=timestamps)
                    modBasesSplit = np.split(modBases,parsedDict['rowBreakIndices'])

#                    def f(fig,axes):
#                        axes.plot(timestamps)
#                    pop(plotFunc=f)
#                    def f(fig,axes):
#                        axes.plot(timestamps,baselines,'rx')
#                        axes.plot(timestamps,np.bitwise_and(modBases,pulseMask),'b.')
#                    pop(plotFunc=f)

                    for iRow in range(len(dataset)):
                        if len(dataset[iRow]) > 0:
                            dataset[iRow] = np.bitwise_and(dataset[iRow],deleteBaseMask) + np.left_shift(np.bitwise_and(modBasesSplit[iRow],pulseMask),nBitsAfterBaseline)
                        

#                    resp = utils.confirm('Continue?')
#                    if resp == False:
#                        exit(0)
    obs.close()
                    

def main():
    run='PAL2012'
    sunsetDate='20121205'
    ts='20121206-035810'
    smoothObs(run,sunsetDate,ts)

if __name__=='__main__':
    main()





