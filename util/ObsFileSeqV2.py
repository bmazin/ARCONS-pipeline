import os
import math
import time
import warnings
import numpy as np
from util import FileName
from util import ObsFile
from util import TCS
from interval import interval
import pyfits
import matplotlib.pyplot as plt
import pickle

from headers.DisplayStackHeaders import writeImageStack, readImageStack


class ObsFileSeq():
    """
    Deal with a sequence of obsFiles, and present data as a set of
    frames.  Each frame has the telescope pointing at the same location,
    as calculated by the TCS log file written by the data acquisition
    system.

    For a given run and date, consider the list of timeStamps as one
    continuous observation, with the given name.  Divide the observations
    into a set of frames, with the maximum amount of seconds in one frame given
    by dt

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

    def __init__(self, name, run, date, timeStamps, dt):
        """
        name -- a useful name for this set of objects
        run -- the data campaign name (LICK2014)
        date -- sundown date (20141020)
        timeStamps -- the UTC date/time stamp
        dt -- the maximum number of seconds for one frame
        beamMapFile -- the beamMap file to use (None to use default)
        """
        self.name = name
        self.run = run
        self.date = date
        self.timeStamps = timeStamps
        self.timeStamps.sort()
        self.dt = dt
        self.obsFiles = []
        self.fileNames = []
        self.obsFileUnixTimes = []
        for timeStamp in self.timeStamps:
            fn = FileName.FileName(run, date, timeStamp)
            self.fileNames.append(fn)
            of = ObsFile.ObsFile(fn.obs(), repeatable=True)
            of.loadBeammapFile(fn.beammap())
            of.loadBestWvlCalFile()
            fn2 = FileName.FileName(run, date, "")
            of.loadFlatCalFile(fn2.flatSoln())
            try:
                of.loadHotPixCalFile(fn.timeMask())
                self.hotPixelsApplied = True
            except:
                self.hotPixelsApplied = False
            self.obsFiles.append(of)
            self.obsFileUnixTimes.append(of.getFromHeader('unixtime'))
        self.tcs = TCS.TCS(run, date)
        self.tcsDict = self.tcs.select(self.obsFiles[0], self.obsFiles[-1])

        # each interval covers one obsFile
        self.obsIntervals = []
        for i in range(len(self.obsFiles)):
            tStart = self.obsFiles[i].getFromHeader('unixtime')
            tEndThis = tStart + self.obsFiles[i].getFromHeader('exptime')
            if i < len(self.obsFiles) - 1:
                tStartNext = self.obsFiles[i + 1].getFromHeader('unixtime')
                tEnd = min(tEndThis, tStartNext)
            else:
                tEnd = tEndThis
            self.obsIntervals.append(interval[tStart, tEnd])
        self._defineFrames(dt)

        # Default settings for astrometry
        try:
            self.setRm()
        except:
            # Fails if the telescope doesn't move every frame...
            # Ignore this for now
            pass

    def setRaDrift(self, raDrift):
        self.raArcsecPerSec = raDrift

    def setScaleThetaDrift(self, moveList, driftList):
        """
        Set scale and theta for this sequence.
    
        input:  
            driftList is a list of dictionaries, containing iFrame, row, col
            for the object found in two frames separated significantly in time
            and with the same ra,dec offset
        
            moveList is a list of dictionaries, containing iFrame, row, col
            for the object found in two frames that are close in time

    """
        # Calculate the scale and rotation with no drift
        matchList = []
        for frc in moveList:
            iFrame = frc['iFrame']
            matchList.append({'ra': self.tcsDict['raOffset'][iFrame] / 3600.0,
                              'dec': self.tcsDict['decOffset'][iFrame] / 3600.0,
                              'row': frc['row'],
                              'col': frc['col']})
        scaleTheta = ObsFileSeq.getScaleTheta(matchList)
        ct = math.cos(scaleTheta['theta'])
        st = math.sin(scaleTheta['theta'])

        # See how much row,col=0,0 moved between the two drift frames
        ras = []
        decs = []
        times = []
        for i in range(2):
            row = driftList[i]['row']
            col = driftList[i]['col']
            print "i,row,col", i, row, col
            ras.append((col * ct + row * st) / scaleTheta['scale'])
            decs.append((-col * st + row * ct) / scaleTheta['scale'])
            iFrame = driftList[i]['iFrame']
            times.append(self.getTimeBySeq(iFrame))
        self.raArcsecPerSec = 3600.0 * (ras[1] - ras[0]) / (times[1] - times[0])
        print "ras=", ras
        print "times=", times
        print "raArcsecPerSec", self.raArcsecPerSec
        self.rc0 = np.zeros((2, len(self.frameObsInfos)), dtype=np.float)
        # Calculate the scale and rotation, including drift
        matchList = []
        t0 = self.getTimeBySeq(0)
        for frc in moveList:
            iFrame = frc['iFrame']
            t = self.getTimeBySeq(iFrame)
            raDrift = self.raArcsecPerSec * (t - t0)
            ra = (self.tcsDict['raOffset'][iFrame] - raDrift) / 3600.0
            matchList.append({'ra': ra,
                              'dec': self.tcsDict['decOffset'][iFrame] / 3600.0,
                              'row': frc['row'],
                              'col': frc['col']})
        scaleTheta = ObsFileSeq.getScaleTheta(matchList)

        # Find the row,col at the ra,dec of each frame
        ct = math.cos(scaleTheta['theta'])
        st = math.sin(scaleTheta['theta'])
        for iFrame in range(len(self.frameObsInfos)):
            t = self.getTimeBySeq(iFrame)
            raDrift = self.raArcsecPerSec * (t - t0)
            print "iFrame,raDrift", iFrame, raDrift
            raOff = (self.tcsDict['raOffset'][iFrame] - raDrift) / 3600.0
            deOff = (self.tcsDict['decOffset'][iFrame]) / 3600.0
            # 0 for row; 1 for col
            self.rc0[0:iFrame] = (raOff * st + deOff * ct) / scaleTheta['scale']
            self.rc0[1:iFrame] = (raOff * ct - deOff * st) / scaleTheta['scale']
            print "iFrame, raOffset, deOffset", iFrame, self.tcsDict['raOffset'][iFrame], self.tcsDict['decOffset'][
                iFrame], raOff, deOff
            # print "iFrame, raOff, deOff, row, col",iFrame,raOff,deOff, self.rc0[0,iFrame], self.rc0[1,iFrame]
        # Numpy Kung-Fu here to subtract the minimum row,col
        self.rc0 -= self.rc0.min(axis=1)[:, None]
        # Calculate the size of the full image to include all pixels
        self.nRowCol = self.rc0.max(axis=1)
        self.nRowCol[0] += self.obsFiles[0].nRow
        self.nRowCol[1] += self.obsFiles[0].nCol
        self.nRowCol = np.ceil(self.nRowCol).astype(np.int)

    @staticmethod
    def getScaleTheta(matchList, flip=1):
        """
        Calculate scale and theta for the measurement of an object in two 
        frames.  The coordinate system is centered on the object.

        col,row = position in pixel coordinates
        ra,dec  = position in sky coordinates

        The transformation is specified with:
        f = +1 or -1 to flip one axis:
        col0,row0 is the location in pixel coordinates of the origin of
            the sky coordinates, where ra,dec = 0,0
        theta -- rotation angle
        scale -- degrees/pixel
        
        The transformation equations are:
        col = col0 + (flip*ra*cos(theta) - dec*sin(theta)) / scale
        row = row0 + (flip*ra*sin(theta) + dec*sin(theta)) / scale
        
        ra  = ( col*cos(theta) + row*sin(theta)) * scale / flip
        dec = (-col*sin(theta) + row*cos(theta)) * scale

        input:
            matchList is a list of dictionaries, containing ra, dec, row, col.
            ra,dec is the location (in decimal degrees) of row,col=(0,0)
            row,col is the location of the object in the frame
            
        return:
            a dictionary of scale (in degrees/pixel) and theta (radians)
        """
        m0 = matchList[0]
        m1 = matchList[1]
        dra = m1['ra'] - m0['ra']
        ddec = m1['dec'] - m0['dec']
        dr = m1['row'] - m0['row']
        dc = m1['col'] - m0['col']
        theta = math.atan2((flip * dra * dr - ddec * dc), (ddec * dr + flip * dra * dc))
        scale = math.sqrt((dra ** 2 + ddec ** 2) / (dc ** 2 + dr ** 2))
        return {"scale": scale, "theta": theta}

    def setTransform(self, scale, thetaDegrees, flip, rdot, cdot):
        """
        """
        self.scale = scale  # degrees/pixel
        self.thetaDegrees = thetaDegrees
        print "in setTransform:  scale, thetaDegrees =",scale, thetaDegrees
        self.ct = math.cos(math.radians(thetaDegrees))
        self.st = math.sin(math.radians(thetaDegrees))
        self.flip = flip
        self.rdot = rdot
        self.cdot = cdot
        # Keep track of the position of r,c=0,0 in each fram with the arreay self.rc0,where
        #r0 = int(self.rc0[0, iFrame])
        #c0 = int(self.rc0[1, iFrame])
        self.rc0 = np.zeros((2, self.nFrames), dtype=float)
        for iFrame in range(self.nFrames):
            r0c0 = self.getR0C0(iFrame)
            self.rc0[0, iFrame] = r0c0['row0']
            self.rc0[1, iFrame] = r0c0['col0']
        # Numpy yoga to subtract minimum row,col and set nRowNcol
        self.rc0 -= self.rc0.min(axis=1)[:, None]
        self.nRowCol = self.rc0.max(axis=1)
        self.nRowCol[0] += self.obsFiles[0].nRow
        self.nRowCol[1] += self.obsFiles[0].nCol
        self.nRowCol = np.ceil(self.nRowCol).astype(np.int)
        print "end of setTransform:  nRowCol=",self.nRowCol

    def getR0C0(self, iFrame):
        ra = self.tcsDict['raOffset'][iFrame]
        dec = self.tcsDict['decOffset'][iFrame]
        col = (self.flip * ra * self.ct - dec * self.st) / self.scale
        row = (self.flip * ra * self.st + dec * self.ct) / self.scale
        dt = self.tcsDict['timeOffset'][iFrame]
        col0 = -col - dt*self.cdot
        row0 = -row - dt*self.rdot
        retval = dict(col0=col0, row0=row0)
        print "iFrame=%3d ra(as)=%5.1f dec(as)=%5.1f row=%5.1f col=%5.1f dt=%9.3f col0=%5.1f row0=%5.1f"%(iFrame, ra, dec, row, col, dt, col0, row0)
        return retval

    def setRm(self,
              degreesPerPixel=0.4 / 3600,
              thetaDeg=0.0,
              raArcsecPerSec=0.0,
              verbose=False):
        """
        Sets variables that will be  used to offset frames
        self.rdl is a list of raOffset, decOffset based on where the 
            telescope says it was pointing, adding in the drift in ra
        self.rc0 is a list of the row,col locations 
        """
        if verbose:
            print " arcsecPerPixel = ", degreesPerPixel * 3600
            print "theta (degrees) = ", thetaDeg
            print " raArcsecPerSec = ", raArcsecPerSec
        self.degreesPerPixel = degreesPerPixel
        self.thetaDeg = thetaDeg
        self.raArcsecPerSec = raArcsecPerSec
        theta = math.radians(thetaDeg)
        sct = math.cos(theta) * degreesPerPixel
        sst = math.sin(theta) * degreesPerPixel
        self.rmPixToEq = np.array([[sct, -sst], [sst, sct]])
        self.rmEqToPix = np.linalg.inv(self.rmPixToEq)
        t0 = self.getTimeBySeq(0)
        self.rdl = []
        for iFrame in range(len(self.frameObsInfos)):
            t = self.getTimeBySeq(iFrame)
            raDrift = raArcsecPerSec * (t - t0)
            raOff = (self.tcsDict['raOffset'][iFrame] - raDrift) / 3600.0
            deOff = (self.tcsDict['decOffset'][iFrame]) / 3600.0
            self.rdl.append([raOff, deOff])
        self.rc0 = np.dot(self.rmEqToPix, np.array(self.rdl).transpose())
        # Numpy Kung-Fu here to subtract the minimum row,col
        self.rc0 -= self.rc0.min(axis=1)[:, None]
        self.nRowCol = self.rc0.max(axis=1)
        self.nRowCol[0] += self.obsFiles[0].nRow
        self.nRowCol[1] += self.obsFiles[0].nCol
        self.nRowCol = np.ceil(self.nRowCol).astype(np.int)

    def makeMosaicImage(self, iFrameList=None, wvlBinRange=None,
                        verbose=False):
        """
        create a mosaic image of the frames listed, in the wavelength bin range

        input:  iFrameList, default None uses all frames
                wvlBinRange, default None uses all wavelength bins,
                otherwise (wBinMin,wBinMax)

        output:  a numpy 2d of the counts/second image
        """
        try:
            self.cubes
        except AttributeError:
            if verbose:
                print "ObsFileSeq.makeMosaicImage:  loadSpectralCubes()"
            self.loadSpectralCubes()
        cubeSum = np.zeros((self.nRowCol[0], self.nRowCol[1]))
        effIntTimeSum = np.zeros((self.nRowCol[0], self.nRowCol[1]))
        nRowCube = self.obsFiles[0].nRow
        nColCube = self.obsFiles[0].nCol
        if iFrameList is None:
            iFrameList = range(self.nFrames)
        if wvlBinRange is None:
            wBinMin = 0
            wBinMax = self.cubes[0]['cube'].shape[2]
        else:
            wBinMin = wvlBinRange[0]
            wBinMax = wvlBinRange[1]
        for iFrame in iFrameList:
            r0 = int(self.rc0[0, iFrame])
            c0 = int(self.rc0[1, iFrame])
            if verbose:
                print "ObsFileSeq:makeMosaicImage:  r0,c0=", r0, c0
            # The third index here is where you select which wavelength bins
            # to include
            cubeSum[r0:r0 + nRowCube, c0:c0 + nColCube] += \
                self.cubes[iFrame]['cube'][:, :, wBinMin:wBinMax].sum(axis=2)
            effIntTimeSum[r0:r0 + nRowCube, c0:c0 + nColCube] += \
                self.cubes[iFrame]['effIntTime'][:, :]
        with np.errstate(divide='ignore'):
            cps = cubeSum / effIntTimeSum
        cps = np.nan_to_num(cps)
        return cps

    def __del__(self):
        for of in self.obsFiles:
            del of

    def _defineFrames(self, dt):
        self.dt = dt
        mts = self.tcsDict['time']
        # make a list of times:
        # start of first observation, each move, end of last observation
        times = np.zeros(len(mts) + 2)
        times[1:-1] = mts
        times[0] = self.obsFiles[0].getFromHeader('unixtime')
        self.beginTime = times[0]
        times[-1] = self.obsFiles[-1].getFromHeader('unixtime') + \
                    self.obsFiles[-1].getFromHeader('exptime')

        # Divide these segments into lengths of dt.  Exclude two seconds after
        # each boundary, to let the telescope settle down after a move
        self.frameIntervals = []
        self.locationIdx = []
        for i in range(len(times) - 1):
            t0 = times[i] + 2
            t1 = times[i + 1]
            nSeg = int((t1 - t0) / dt) + 1
            delta = (t1 - t0) / nSeg
            for j in range(nSeg):
                self.frameIntervals.append(interval[t0 + j * delta,
                                                    t0 + (j + 1) * delta])
                self.locationIdx.append(i)
        # For each frame, determine a list of:
        #      obsFile, firstSec, integrationTime
        # In general, one frame will need more than one of these if it
        # spans the time boundary between obsFiles
        self.frameObsInfos = []
        for iFrame in range(len(self.frameIntervals)):
            frameObsInfo = []
            thisInterval = self.frameIntervals[iFrame]
            # see if this interval overlaps with an obsFile
            for i, obsInterval in enumerate(self.obsIntervals):
                overlap = thisInterval & obsInterval
                if len(overlap) > 0:
                    tBeg = overlap[0][0]
                    tEnd = overlap[0][1]
                    integrationTime = tEnd - tBeg
                    firstSec = tBeg - \
                               self.obsFiles[i].getFromHeader('unixtime')
                    obs = self.obsFiles[i]
                    obsInfo = {"obs": obs,
                               "iObs": i,
                               "firstSec": firstSec,
                               "integrationTime": integrationTime}
                    frameObsInfo.append(obsInfo)
            self.frameObsInfos.append(frameObsInfo)
        self.nFrames = len(self.frameObsInfos)

    def getTimeBySeq(self, iSeq):
        """
        get the mean time of the frame
        """
        foi = self.frameObsInfos[iSeq]
        wtSum = 0
        wSum = 0
        for oi in foi:
            w = oi['integrationTime']
            t = self.obsFileUnixTimes[oi['iObs']] + oi['firstSec'] + \
                0.5 * oi['integrationTime']
            wtSum += w * t
            wSum += w
        meanTime = wtSum / float(wSum)
        return meanTime

    def getTargetList(self, printLines=True):
        """
        get list of information:  the timstamp + target description in header
        default printLine=True to also print each line to stdout

        return:  a list of the lines
        """
        retval = []
        for i, timeStamp in enumerate(self.timeStamps):
            target = self.obsFiles[i].getFromHeader('target')
            line = "%2d %s %s" % (i, timeStamp, target)
            if printLines:
                print line
            retval.append(line)
        return retval

    def getFrameList(self, printLines=True):
        """
        returns a list of lines which describes the frames.
        Each line has an index, the time (relative to the time of the first
        frame), effective integration time, location number,
        the (raOffset,decOffset), and a list of obs files used
        in the frame.
        """
        retval = []
        for i, frameInterval in enumerate(self.frameIntervals):
            t0 = frameInterval[0][0] - self.beginTime
            t1 = frameInterval[0][1] - self.beginTime
            dt = t1 - t0
            locIdx = self.locationIdx[i]
            xOff = self.tcsDict['raOffset'][locIdx]
            yOff = self.tcsDict['decOffset'][locIdx]
            obsFiles = ""
            for frameObsInfo in self.frameObsInfos[i]:
                obsFiles += " %d" % frameObsInfo['iObs']
            fmt = "i=%3d  begin=%8.2f expTime=%6.2f loc=%2d  (%5.1f,%5.1f) %s"
            line = fmt % (i, t0, dt, locIdx, xOff, yOff, obsFiles)
            if printLines:
                print line
            retval.append(line)
        return retval

    def getFrameDict(self, printLines=True):
        """
        Pretty much the same as getFramList except that it returns
        a dictionary for each frame
        """
        frameInfo = []
        for i, frameInterval in enumerate(self.frameIntervals):
            t0 = frameInterval[0][0] - self.beginTime
            t1 = frameInterval[0][1] - self.beginTime
            dt = t1 - t0
            locIdx = self.locationIdx[i]
            xOff = self.tcsDict['raOffset'][locIdx]
            yOff = self.tcsDict['decOffset'][locIdx]
            meanTF = self.getTimeBySeq(i)
            fI = {"iframe": i,
                  "begin": t0,
                  "expTime": dt,
                  "loc": locIdx,
                  "offsRA": xOff,
                  "offsDec": yOff,
                  "meanTime": meanTF}
            # "obsFile":ofs.fileNames[i].obs(),
            # "ob":ofs.obsFiles[i]}
            if printLines:
                print fI
            frameInfo.append(fI)

        self.frameDict = frameInfo
        return frameInfo

    def getSpectralCubeByFrame(self, iFrame, weighted=False,
                               fluxWeighted=False,
                               wvlStart=None, wvlStop=None,
                               wvlBinWidth=None, energyBinWidth=None,
                               wvlBinEdges=None, timeSpacingCut=None):
        """
        return the spectral cube for this frame

        call ObsFile.getSpectralCube for each ObsFile in this frame.
        The dictionary returned copies 'wvlBinEdges' from the first ObsFile,
        and sums the 'cube' and 'effIntTime' from all ObsFiles.

        I left the print statements in to report progress, because this is
        very slow.

        """
        retval = None
        thisInterval = self.frameIntervals[iFrame]
        for i, ofInterval in enumerate(self.obsIntervals):
            overlap = thisInterval & ofInterval
            if len(overlap) > 0:
                tBeg = overlap[0][0]
                tEnd = overlap[0][1]
                integrationTime = tEnd - tBeg
                firstSec = tBeg - self.obsFiles[i].getFromHeader('unixtime')
                obs = self.obsFiles[i]
                obs.setWvlCutoffs(wvlLowerLimit=wvlStart,
                                  wvlUpperLimit=wvlStop)
                spectralCube = \
                    obs.getSpectralCube(firstSec=firstSec,
                                        integrationTime=integrationTime,
                                        weighted=weighted,
                                        fluxWeighted=fluxWeighted,
                                        wvlStart=wvlStart,
                                        wvlStop=wvlStop,
                                        wvlBinWidth=wvlBinWidth,
                                        energyBinWidth=energyBinWidth,
                                        wvlBinEdges=wvlBinEdges,
                                        timeSpacingCut=timeSpacingCut
                                        )
                cube = spectralCube['cube']
                wbe = spectralCube['wvlBinEdges']
                eit = spectralCube['effIntTime']
                if retval is None:
                    retval = {'cube': cube,
                              'wvlBinEdges': wbe,
                              'effIntTime': eit}
                else:
                    retval['cube'] += cube
                    retval['effIntTime'] += eit
        return retval

    def loadSpectralCubes(self, weighted=False, fluxWeighted=False,
                          wvlStart=None, wvlStop=None,
                          wvlBinWidth=None, energyBinWidth=None,
                          wvlBinEdges=None, timeSpacingCut=None):
        """
        calls getSpectralCubeByFrame on each iFrame, storing the
        results in the list self.cubes

        use a pickle file named name.pkl as a buffer.  If that file
        exists, load the cubes from there, and save the cubes there
        after loading.  WARNING -- the settings are not stored, so
        they are ignored when loading from the pickle file.

        """
        cpfn = self.name + "-cubes.pkl"
        if os.path.isfile(cpfn):
            print "loadSpectralCubes:  load from ", cpfn
            self.cubes = pickle.load(open(cpfn, 'rb'))
        else:
            self.cubes = []
            for iFrame in range(len(self.frameIntervals)):
                print "now load spectral cube for iFrame=", iFrame
                cube = self.getSpectralCubeByFrame(iFrame)
                print "counts are ", cube['cube'].sum()
                self.cubes.append(cube)
                #                self.cubes.append(self.getSpectralCubeByFrame(iFrame,
                #                                                              weighted,
                #                                                              fluxWeighted,
                #                                                             wvlStart,
                #                                                             wvlStop,
                #                                                             wvlBinWidth,
                #                                                             energyBinWidth,
                #                                                             wvlBinEdges,
                #                                                             timeSpacingCut))
                print "counts read:  ", self.cubes[-1]['cube'].sum()
            pickle.dump(self.cubes, open(cpfn, 'wb'))

    def makePngFileByInterval(self, thisInterval, wvMin=3000, wvMax=12000,
                              rateMax=None):
        fn = "%s-%03d-%05d-%05d.png" % \
             (self.name, thisInterval, int(wvMin), int(wvMax))
        print "now make fn=", fn
        cubeSum = self.cubes[thisInterval]['cube'].sum(axis=2)
        effIntTimeSum = self.cubes[thisInterval]['effIntTime']
        old_settings = np.seterr(all='ignore')
        np.seterr(divide='ignore')
        rate = np.nan_to_num(cubeSum / effIntTimeSum)
        np.seterr(**old_settings)
        print fn, rate.min(), rate.max()
        plt.clf()
        if rateMax is None:
            rateMax = rate.max()
        plt.pcolor(rate, cmap='hot', vmin=0, vmax=rateMax)
        plt.colorbar()
        try:
            os.remove(fn)
        except OSError:
            pass
        plt.title(fn)
        plt.savefig(fn)

    def makeAllFitsFiles(self, wvMin, wvMax):
        self.loadSpectralCubes(wvlStart=wvMin, wvlStop=wvMax)
        for interval in range(len(self.frameIntervals)):
            self.makeFitsFileByInterval(interval, wvMin, wvMax)
        return

    def makeFitsFileByInterval(self, thisInterval, wvMin, wvMax):
        fn = "%s-%03d-%05d-%05d.fit" % (self.name, thisInterval,
                                        int(wvMin), int(wvMax))
        print "now make fn=", fn
        pixels = self.cubes[thisInterval]['cube'].sum(axis=2)
        print "number of counts=", pixels.sum()

        hdu = pyfits.PrimaryHDU(pixels)
        try:
            os.remove(fn)
        except OSError:
            pass
        hdu.writeto(fn)

    def loadImageStack(self, fileName, wvlStart=None, wvlStop=None,
                       weighted=True, fluxWeighted=False,
                       getRawCount=False, scaleByEffInt=True,
                       deadTime=100.e-6):
        # If the file exists, read it out
        if os.path.isfile(fileName):
            return readImageStack(fileName)
        # if the file doesn't exists, make it
        else:
            images = []
            pixIntTimes = []
            startTimes = []
            endTimes = []
            intTimes = []
            for iFrame in range(len(self.frameIntervals)):
                im_dict = self.getPixelCountImageByFrame(iFrame,
                                                         wvlStart, wvlStop,
                                                         weighted,
                                                         fluxWeighted,
                                                         getRawCount,
                                                         scaleByEffInt,
                                                         deadTime)
                images.append(im_dict['image'])
                pixIntTimes.append(im_dict['pixIntTime'])
                startTimes.append(im_dict['startTime'])
                endTimes.append(im_dict['endTime'])
                intTimes.append(im_dict['intTime'])

            writeImageStack(fileName, images, startTimes=startTimes,
                            endTimes=endTimes, intTimes=intTimes,
                            pixIntTimes=pixIntTimes, targetName=self.name,
                            run=self.run,
                            nFrames=len(self.frameIntervals),
                            wvlLowerLimit=wvlStart,
                            wvlUpperLimit=wvlStop, weighted=weighted,
                            fluxWeighted=fluxWeighted,
                            hotPixelsApplied=self.hotPixelsApplied,
                            maxExposureTime=self.dt,
                            tStamps=self.timeStamps)
            # return {'images':images,'pixIntTimes':pixIntTimes,
            # 'startTimes':startTimes,'endTimes':endTimes,'intTimes':intTimes}
            return readImageStack(fileName)

    def getPixelCountImageByFrame(self, iFrame, wvlStart=None, wvlStop=None,
                                  weighted=True, fluxWeighted=True,
                                  getRawCount=False, scaleByEffInt=True,
                                  deadTime=100.e-6):
        '''
        This gets the i'th image

        Inputs:
            iFrame - which frame you want
            wvlStart - starting wavelength range
            wvlStop - ending wavelength range
            weighted, fluxWeighted, getRawCount, scaleByEffInt -
                options for obsFile.getPixelCountImage()
            deadTime - for deadtime correcting image

        Returns:
            Dictionary with the following keys:
            'image' - fully calibrated, corrected image.
                      scaled to the total integration time
                      deadTime corrected
            'pixIntTime' - actual integration time for each pixel in image
            'intTime' - length of exposure
            'startTime' - beginning of image (unix time)
            'endTime' - end of image. Might be different from
                startTime+intTime if there's a break in the middle while
                switching to a new obsFile

        '''

        retval = None
        for obsInfo in self.frameObsInfos[iFrame]:
            print obsInfo
            obsInfo['obs'].setWvlCutoffs(wvlLowerLimit=wvlStart,
                                         wvlUpperLimit=wvlStop)
            im_dict = obsInfo['obs']. \
                getPixelCountImage(obsInfo["firstSec"],
                                   obsInfo["integrationTime"],
                                   weighted,
                                   fluxWeighted,
                                   getRawCount,
                                   scaleByEffInt=False)
            # Do this manually so
            # we can correct deadTime first
            im = im_dict['image']
            # print 'im: ', np.sum(im)
            # Correct for dead time
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",
                                        'invalid value encountered in divide',
                                        RuntimeWarning)
                w_deadTime = \
                    1.0 - im_dict['rawCounts'] * deadTime / im_dict['effIntTimes']
            im = im / w_deadTime
            if scaleByEffInt:
                # Correct for exposure time
                im = im * obsInfo["integrationTime"] / im_dict['effIntTimes']
            # Remove any funny values
            im[np.invert(np.isfinite(im))] = 0.
            # print '--> ', np.sum(im)
            if retval is None:
                retval = {'image': im, 'pixIntTime': im_dict['effIntTimes'],
                          'intTime': obsInfo["integrationTime"],
                          'startTime': self.frameIntervals[iFrame][0][0],
                          'endTime': self.frameIntervals[iFrame][0][1]}
            else:
                retval['image'] += im
                retval['pixIntTime'] += im_dict['effIntTimes']
                retval['intTime'] += obsInfo["integrationTime"]
        return retval

    def plotLocations(self, fileName=None):
        plt.clf()
        x = self.tcsDict['raOffset']
        y = self.tcsDict['decOffset']
        plt.plot(x, y)
        for i in range(len(x)):
            plt.text(x[i], y[i], str(i),
                     horizontalalignment="center",
                     verticalalignment="center")
        plt.axes().set_aspect('equal', 'datalim')
        plt.title(self.name)
        plt.xlabel("raOffset (arcsec)")
        plt.ylabel("decOffset (arcsec)")
        print "in ObsFileSeq.plotLocations:  fileName=", fileName
        if not fileName:
            plt.show()
        else:
            plt.savefig(fileName)


if __name__ == "__main__":
    if 0:
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
        ofs = ObsFileSeq(name, run, date, tsl, dt)
        print "Now call getTargetList"
        ofs.getTargetList()
        print "Now call getFrameList"
        ofs.getFrameList()
        ofs.plotLocations(name + ".png")
        print "now get time of first frame"
        for i in range(66):
            print "i=", i, " time=", ofs.getTimeBySeq(i)
        # apci = ofs.getAllPixelCountImages(getRawCount=True)
        del ofs
