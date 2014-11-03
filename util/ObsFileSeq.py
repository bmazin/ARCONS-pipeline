import os
import numpy as np
from util import FileName
from util import ObsFile
from util import TCS
from interval import interval, inf, imath
import pyfits
import matplotlib.pyplot as plt
class ObsFileSeq():
    """
    Deal with a sequence of obsFiles

    For a given run and date, consider the list of timeStamps as one
    continuous observation, with the given name.  Divide the observations
    into a set of frames, with the maximum amount of seconds in one frame given
    by dt.

    If the data acquisition system was doing a mosaic at the time, each
    frame will contain data only for one telescope position.  Two seconds after
    each telescope move is masked, to give the telescope a chance to settle 
    in its new position.

    One frame may contain data from more than one ObsFile.

    getTargetList lists the timestamp (and target description) to show
    the list of obsFiles included in this sequence.

    getFrameList lists the information for each frame;  its ra,dec offset,
    and which obsFiles contribute to this frame.

    plotLocations makes a simple png of the frame numbers plotted in 
    raOffset,decOffset from the initial position.

    executing just this file demonstrates calls to getTargetlist and
    plotLocations for a mosaic of the ring nebula.


    """
    def __init__(self,name,run,date,timeStamps,dt):
        """
        name -- a useful name for this set of objects
        run -- the data campaign name (LICK2014)
        date -- sundown date (20141020)
        timeStamps -- the UTC date/time stamp
        dt -- the maximum number of seconds for one frame
        """
        self.name = name
        self.run = run
        self.date = date
        self.timeStamps = timeStamps
        self.timeStamps.sort()
        self.obsFiles = []
        self.fileNames = []
        for timeStamp in self.timeStamps:
            fn = FileName.FileName(run,date,timeStamp)
            self.fileNames.append(fn)
            of = ObsFile.ObsFile(fn.obs()) 
            self.obsFiles.append(of)
        self.tcs = TCS.TCS(run,date)
        self.tcsDict = self.tcs.select(self.obsFiles[0],self.obsFiles[-1])

        # each interval covers one obsFile
        self.obsIntervals = []
        for i in range(len(self.obsFiles)):
            tStart = self.obsFiles[i].getFromHeader('unixtime')
            tEndThis = tStart+self.obsFiles[i].getFromHeader('exptime')
            if i < len(self.obsFiles)-1:
                tStartNext = self.obsFiles[i+1].getFromHeader('unixtime')
                tEnd = min(tEndThis, tStartNext)
            else:
                tEnd = tEndThis
            self.obsIntervals.append(interval[tStart,tEnd])
        self._defineFrames(dt)
    def __del__(self):
        for of in self.obsFiles:
            del of
    def _defineFrames(self,dt):
        self.dt = dt
        mts = self.tcsDict['time']
        print "number of moves=",len(mts)
        # make a list of times: 
        #start of first observation, each move, end of last observation
        times = np.zeros(len(mts)+2) 
        times[1:-1] = mts
        times[0]  = self.obsFiles[0].getFromHeader('unixtime')        
        self.beginTime = times[0]
        times[-1] = self.obsFiles[-1].getFromHeader('unixtime') + self.obsFiles[-1].getFromHeader('exptime') 

        # Divide these segments into lengths of dt.  Exclude two seconds after
        # each boundary, to let the telescope settle down after a move
        self.frameIntervals = []
        self.locationIdx = []
        for i in range(len(times)-1):
            t0 = times[i]+2
            t1 = times[i+1]
            nSeg = int((t1-t0)/dt)+1
            delta = (t1-t0)/nSeg
            for j in range(nSeg):
                self.frameIntervals.append(interval[t0+j*delta,
                                                    t0+(j+1)*delta])
                self.locationIdx.append(i) 
        # For each frame, determine a list of:  obsFile, firstSec, integrationTime
        # In general, one frame will need more than one of these if it
        # spans the time boundary between obsFiles
        self.frameObsInfos = []
        for iFrame in range(len(self.frameIntervals)):
            frameObsInfo = []
            thisInterval = self.frameIntervals[iFrame]
            # see if this interval overlaps with an obsFile
            for i,obsInterval in enumerate(self.obsIntervals):
                overlap = thisInterval & obsInterval
                if len(overlap) > 0:
                    tBeg = overlap[0][0]
                    tEnd = overlap[0][1]
                    integrationTime = tEnd-tBeg
                    firstSec = tBeg - self.obsFiles[i].getFromHeader('unixtime')
                    obs = self.obsFiles[i]
                    obsInfo = {"obs":obs,
                               "iObs":i,
                               "firstSec":firstSec,
                               "integrationTime":integrationTime}
                    frameObsInfo.append(obsInfo)
            self.frameObsInfos.append(frameObsInfo)

    def getTargetList(self,printLines=True):
        """
        get list of information:  the timstamp + target description in header
        default printLine=True to also print each line to stdout
        
        return:  a list of the lines
        """
        retval = []
        for i,timeStamp in enumerate(self.timeStamps):
            target = self.obsFiles[i].getFromHeader('target')
            line = "%2d %s %s"%(i,timeStamp,target)
            if printLines:
                print line
            retval.append(line)
        return retval
    def getFrameList(self,printLines=True):
        """
        returns a list of lines which describes the frames.
        Each line has an index, the time (relative to the time of the first
        frame), effective integration time, location number,
        the (raOffset,decOffset), and a list of obs files used
        in the frame.  
        """
        retval = []
        for i,frameInterval in enumerate(self.frameIntervals):
            t0 = frameInterval[0][0]-self.beginTime
            t1 = frameInterval[0][1]-self.beginTime
            dt = t1-t0
            locIdx = self.locationIdx[i]
            xOff = self.tcsDict['raOffset'][locIdx]
            yOff = self.tcsDict['decOffset'][locIdx]
            obsFiles = ""
            for frameObsInfo in self.frameObsInfos[i]:
                obsFiles += " %d"%frameObsInfo['iObs']
            line = "i=%3d  begin=%8.2f expTime=%6.2f loc=%2d  (%5.1f,%5.1f) %s"%(i,t0,dt,locIdx,xOff,yOff,obsFiles)
            if printLines:
                print line
            retval.append(line)
        return retval


    def getSpectralCubes(self,iFrame,wvMin,wvMax):
        thisInterval = self.frameIntervals[iFrame]
        return self.getSpectralCubesByInterval(thisInterval,wvMin,wvMax)

    def getSpectralCubesByInterval(self,thisInterval,wvMin,wvMax):
        retval = []
        for i,ofInterval in enumerate(self.obsIntervals):
            overlap = thisInterval & ofInterval
            if len(overlap) > 0:
                tBeg = overlap[0][0]
                tEnd = overlap[0][1]
                integrationTime = tEnd-tBeg
                firstSec = tBeg - self.obsFiles[i].getFromHeader('unixtime')
                print i,firstSec,integrationTime
                obs = self.obsFiles[i]
                obs.loadBestWvlCalFile() 
                obs.setWvlCutoffs(wvlLowerLimit=wvMin, wvlUpperLimit=wvMax)
                hotPixCalFileName = self.fileNames[i].timeMask()
                obs.loadHotPixCalFile(hotPixCalFileName)
                #obs.loadFlatCalFile(flatCalFileName)
                #obs.loadFluxCalFile(fluxCalFileName)

                # now get the scInfo: dict of cube,wvlBinEdges,effIntTime
                spectralCube = obs.getSpectralCube(firstSec=firstSec,
                                                   integrationTime=integrationTime,
                                                   weighted=False,
                                                   wvlBinWidth=1000
                                                   )
                
                retval.append(spectralCube)
        return retval
    def makePngFileByInterval(self,thisInterval,wvMin,wvMax,maxRate):
        fn = "%s-%03d-%05d-%05d.png"%(self.name,thisInterval,int(wvMin),int(wvMax))
        print "now make fn=",fn
        spectralCubes = self.getSpectralCubes(thisInterval,wvMin,wvMax)
        cubeSum = None
        for spectralCube in spectralCubes:
            cube = spectralCube['cube'].sum(axis=2)
            effIntTime = spectralCube['effIntTime']
            if cubeSum is None:
                cubeSum = cube
                effIntTimeSum = effIntTime
            else:
                cubeSum += cube
                effIntTimeSum += effIntTime
        old_settings = np.seterr(all='ignore')
        np.seterr(divide='ignore')
        rate = np.nan_to_num(cubeSum/effIntTimeSum)
        np.seterr(**old_settings)
        print fn,rate.min(),rate.max()
        plt.clf()
        plt.pcolor(rate, cmap='hot', vmin=0,vmax=maxRate)
        plt.colorbar()
        try:
            os.remove(fn)
        except OSError:
            pass
        plt.title(fn)
        plt.savefig(fn)

    def makeFitsFileByInterval(self,thisInterval,wvMin,wvMax):
        fn = "%s-%03d-%05d-%05d.fit"%(self.name,thisInterval,int(wvMin),int(wvMax))
        print "now make fn=",fn
        scInfo = self.getSpectralCubes(thisInterval,wvMin,wvMax)
        pixels = scInfo[0]['cube'].sum(axis=2)
        print "number of counts=",pixels.sum()
        
        hdu = pyfits.PrimaryHDU(pixels)
        try:
            os.remove(fn)
        except OSError:
            pass
        hdu.writeto(fn)
    def getPixelCountImage(self,t0,t1,weighted=False,
                           fluxWeighted=False, getRawCount=False,
                           scaleByEffInt=False):
        thisInterval = interval[t0,t1]
        for i,ofInterval in enumerate(self.obsIntervals):
            overlap = thisInterval & ofInterval
            if len(overlap) > 0:
                tBeg = overlap[0][0]
                tEnd = overlap[0][1]
                integrationTime = tEnd-tBeg
                firstSec = tBeg - self.obsFiles[i].getFromHeader('unixtime')
                print i,firstSec,integrationTime
                pci = self.obsFiles[i].getPixelCountImage(firstSec,
                                                          integrationTime,
                                                          weighted,
                                                          fluxWeighted,
                                                          getRawCount,
                                                          scaleByEffInt)
        return None

    def plotLocations(self,fileName=None):
        plt.clf()
        x = self.tcsDict['raOffset']
        y = self.tcsDict['decOffset']
        plt.plot(x,y)
        for i in range(len(x)):
            plt.text(x[i],y[i],str(i),
                     horizontalalignment="center",
                     verticalalignment="center")
        plt.axes().set_aspect('equal', 'datalim')
        plt.title(self.name)
        plt.xlabel("raOffset (arcsec)")
        plt.ylabel("decOffset (arcsec)")
        if not fileName:
            plt.show()
        else:
            plt.savefig(fileName)

if __name__ == "__main__":
    name = 'ring-20141020'
    run = "PAL2014"
    date = "20141020"
    tsl = [
        '20141021-033954',
        '20141021-034532',
        '20141021-035035',
        '20141021-035538',
        '20141021-040041',
        '20141021-040544',
        '20141021-041047',
        ]
    dt = 200
    ofs = ObsFileSeq(name,run,date,tsl,dt)
    print "Now call getTargetList"
    ofs.getTargetList()
    print "Now call getFrameList"
    ofs.getFrameList()
    ofs.plotLocations(name+".png")
    #apci = ofs.getAllPixelCountImages(getRawCount=True)
    del ofs
