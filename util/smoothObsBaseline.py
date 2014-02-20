import numpy as np
from util.FileName import FileName
from util import utils
from util.popup import plotArray,pop
import scipy.stats
import tables
import os
import re
import shutil
import sys
from util.readDict import readDict

class BaseSmoother:
    tickDuration = 1.e-6 #s
    nBitsInChannel = 8
    nBitsInParabolaPulseHeight = 12
    nBitsInPulseHeight = 12
    nBitsInBaseline = 12
    nBitsInTimestamp = 20

    nBitsAfterBaseline = 20

    def __init__(self,paramFile):
        #  8 bits - channel
        # 12 bits - Parabola Fit Peak Height
        # 12 bits - Sampled Peak Height
        # 12 bits - Low pass filter baseline
        # 20 bits - Microsecond timestamp
        #bitmask of 12 ones
        self.pulseMask = int(self.nBitsInPulseHeight * '1', 2)
        #bitmask of 20 ones
        self.timestampMask = int(self.nBitsInTimestamp * '1', 2)
        self.deleteBaseMask = int((8+12+12)*'1'+12*'0'+20*'1',2)

        self.params = readDict()
        self.params.read_from_file(paramFile)
        self.run = self.params['run']
        self.smoothRun = self.run+'s'
        self.timeToAverage = self.params['timeToAverage']
        self.trimFraction = self.params['trimFraction']
        self.sunsetDate = self.params['sunsetDate']
        self.obsSequence = self.params['obsSequence']
        self.obsFNs = [FileName(run=self.run,date=self.sunsetDate,tstamp=obsTstamp) for obsTstamp in self.obsSequence]
        self.smoothObsFNs = [FileName(run=self.smoothRun,date=self.sunsetDate,tstamp=obsTstamp) for obsTstamp in self.obsSequence]



    def parsePhotonPackets(self,packets):
            """
            Parses an array of uint64 packets with the obs file format
            inter is an interval of time values to mask out
            returns a list of timestamps,parabolaFitPeaks,baselines
            """
            # parse all packets
            packetsAll = np.array(packets, dtype='uint64') #64 bit photon packet
            timestamps = np.bitwise_and(packets, self.timestampMask)
            baselines = np.bitwise_and(np.right_shift(packets, self.nBitsAfterBaseline), self.pulseMask)
            return {'timestamps':timestamps, 'baselines':baselines}

    def getTimedPacketList(self,pixelData):
        timestamps = []
        baselines = []
        nRows = len(pixelData)
        rowBreaks = np.zeros(nRows-1)
        breakIndex = 0
        for t in range(nRows):
                parseDict = self.parsePhotonPackets(pixelData[t])
                times = parseDict['timestamps']
                bases = parseDict['baselines']
                times = self.tickDuration * times + t

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

    def smoothBaselines(self,baselines,timestamps,timeToAverage=500e-3,trimFraction=0.5):
        #A time based smoothing
        #For each baseline point, all points within a given time interval around the point is
        #trimmed and averaged to make a new baseline point
        modBases = np.array(baselines,dtype=np.uint32)
        for i,ts in enumerate(timestamps):
            startIdx = np.searchsorted(timestamps,ts-timeToAverage/2.)
            endIdx = np.searchsorted(timestamps,ts+timeToAverage/2.)
            trimmedSample = scipy.stats.mstats.trim(baselines[startIdx:endIdx],(trimFraction,0),relative=True)
            if np.all(trimmedSample.mask) == False: #if there are some points that are not trimmed, take the average
                modBases[i] = int(np.ma.mean(trimmedSample))
        return modBases

    def smoothAllObs(self):
        self.obsFileNames = []
        for fn,smoothFN in zip(self.obsFNs,self.smoothObsFNs):#check for both obs and cal files
            if os.path.exists(fn.obs()):
                print 'Copying',fn.obs()
                shutil.copyfile(fn.obs(),smoothFN.obs())
                print 'Smoothing',smoothFN.obs()
                self.smoothObs(smoothFN.obs(),timeToAverage=self.timeToAverage,trimFraction=self.trimFraction)

            elif os.path.exists(fn.cal()):
                print 'Copying',fn.cal()
                shutil.copyfile(fn.cal(),smoothFN.cal())
                print 'Smoothing',smoothFN.cal()
                self.smoothObs(smoothFN.cal(),timeToAverage=self.timeToAverage,trimFraction=self.trimFraction)
            else:
                print 'ERROR: file',fn.obs(),'does not exist, skipping'

    def smoothObs(self,obsFileName,timeToAverage=5.e-3,trimFraction=0.5):

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
                    #print pixelGroup._v_name,
                    for dataset in pixelGroup: #There is only one
                        #print len(dataset),
                        data = dataset.read()
                        parsedDict = self.getTimedPacketList(data)
                        timestamps = parsedDict['timestamps']
                        baselines = parsedDict['baselines']
                        #print len(timestamps)
                        modBases = self.smoothBaselines(baselines=baselines,timestamps=timestamps,timeToAverage=timeToAverage,trimFraction=trimFraction)
                        modBasesSplit = np.split(modBases,parsedDict['rowBreakIndices'])

    #                    def f(fig,axes):
    #                        axes.plot(timestamps)
    #                    pop(plotFunc=f)
    #                    def f(fig,axes):
    #                        axes.plot(timestamps,baselines,'rx-')
    #                        axes.plot(timestamps,np.bitwise_and(modBases,self.pulseMask),'b.-')
    #                    pop(plotFunc=f)
    #
                        for iRow in range(len(dataset)):
                            if len(dataset[iRow]) > 0:
                                dataset[iRow] = np.bitwise_and(dataset[iRow],self.deleteBaseMask) + np.left_shift(np.bitwise_and(modBasesSplit[iRow],self.pulseMask),self.nBitsAfterBaseline)
                            

    #                    resp = utils.confirm('Continue?')
    #                    if resp == False:
    #                        exit(0)
        obs.close()
                    

def main():
    run='PAL2012'
    sunsetDate='20121210'
    ts='20121211-135526'
    print 'starting',ts


if __name__ == '__main__':
    paramFile = sys.argv[1]
    baseSmoother = BaseSmoother(paramFile)
    baseSmoother.smoothAllObs()



