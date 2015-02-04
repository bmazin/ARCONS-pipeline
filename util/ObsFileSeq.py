import os, math, time
import numpy as np
from util import FileName
from util import ObsFile
from util import TCS
from interval import interval, inf, imath
import pyfits
import matplotlib.pyplot as plt
import pickle

from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4.QtGui import *
import astrometry.CentroidCalc as cc
import matplotlib as mpl

class ObsFileSeq():
    """
    Deal with a sequence of obsFiles, and present data as a set of
    frames.  Each frame has the telescope pointing at the same location,
    as calculated by the TCS log file written by the data acquisition
    system.

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
        self.obsFileUnixTimes = []
        for timeStamp in self.timeStamps:
            fn = FileName.FileName(run,date,timeStamp)
            self.fileNames.append(fn)
            of = ObsFile.ObsFile(fn.obs()) 
            fn2 = FileName.FileName(run,date,"")
            of.loadBestWvlCalFile()
            of.loadFlatCalFile(fn2.flatSoln())
            of.loadHotPixCalFile(fn.timeMask())
            self.obsFiles.append(of)
            self.obsFileUnixTimes.append(of.getFromHeader('unixtime'))
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

        # Default settings for astrometry
        self.setRm()
    def setRm(self,
              degreesPerPixel = 0.4/3600,
              thetaDeg = 0.0,
              raArcsecPerSec = 0.0):
        self.degreesPerPixel = degreesPerPixel
        self.thetaDeg = thetaDeg
        self.raArcsecPerSec = raArcsecPerSec
        theta = math.radians(thetaDeg)
        sct = math.cos(theta)*degreesPerPixel
        sst = math.sin(theta)*degreesPerPixel
        self.rmPixToEq = np.array([[sct,-sst],[sst,sct]])
        self.rmEqToPix = np.linalg.inv(self.rmPixToEq)
        t0 = self.getTimeBySeq(0)
        self.rdl = []
        for iFrame in range(len(self.frameObsInfos)):
            t = self.getTimeBySeq(iFrame)
            raDrift = raArcsecPerSec*(t-t0)
            raOff = (self.tcsDict['raOffset'][iFrame]-raDrift)/3600.0
            deOff = (self.tcsDict['decOffset'][iFrame])/3600.0
            self.rdl.append([raOff,deOff])
        self.rc0 = np.dot(self.rmEqToPix,np.array(self.rdl).transpose())
        # Numpy Kung-Fu here to subtract the minimum row,col
        self.rc0 -= self.rc0.min(axis=1)[:,None]
        self.nRowCol = self.rc0.max(axis=1)
        self.nRowCol[0] += self.obsFiles[0].nRow
        self.nRowCol[1] += self.obsFiles[0].nCol
        self.nRowCol = np.ceil(self.nRowCol).astype(np.int)

    def makeMosaicImage(self,iFrameList=None, wvlBinRange=None):
        """
        create a mosaic image of the frames listed, in the wavelength bin range

        input:  iFrameList, default None uses all frames
                wvlBinRange, default None uses all wavelength bins, otherwise (wBinMin,wBinMax)

        output:  a numpy 2d of the counts/second image
        """
        cubeSum = np.zeros((self.nRowCol[0],self.nRowCol[1]))
        effIntTimeSum = np.zeros((self.nRowCol[0],self.nRowCol[1]))
        nRowCube = self.obsFiles[0].nRow
        nColCube = self.obsFiles[0].nCol
        if iFrameList is None:
            iFrameList = range(len(self.frameObsInfos))
        if wvlBinRange is None:
            wBinMin = 0
            wBinMax = self.cubes[0]['cube'].shape[2]
        else:
            wBinMin = wvlBinRange[0]
            wBinMax = wvlBinRange[1]
        for iFrame in iFrameList:
            r0 = int(self.rc0[0,iFrame])
            c0 = int(self.rc0[1,iFrame])
            # The third index here is where you select which wavelength bins to include
            cubeSum[r0:r0+nRowCube,c0:c0+nColCube] += self.cubes[iFrame]['cube'][:,:,wBinMin:wBinMax].sum(axis=2)
            effIntTimeSum[r0:r0+nRowCube,c0:c0+nColCube] += self.cubes[iFrame]['effIntTime'][:,:]
        with np.errstate(divide='ignore'):
            cps = cubeSum/effIntTimeSum
        cps = np.nan_to_num(cps)
        return cps
    def __del__(self):
        for of in self.obsFiles:
            del of

    def _defineFrames(self,dt):
        self.dt = dt
        mts = self.tcsDict['time']
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

    def getTimeBySeq(self,iSeq):
        """
        get the mean time of the frame
        """
        foi = self.frameObsInfos[iSeq]
        wtSum = 0
        wSum = 0
        for oi in foi:
            w = oi['integrationTime']
            t = self.obsFileUnixTimes[oi['iObs']]+oi['firstSec'] + 0.5*oi['integrationTime']
            wtSum += w*t
            wSum += w
        meanTime = wtSum/float(wSum)
        return meanTime

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

    def getFrameDict(self, printLines = True):
        """
        Pretty much the same as getFramList except that it returns a dictionary for each frame
        """
        frameInfo = []
        for i,frameInterval in enumerate(self.frameIntervals):
            t0 = frameInterval[0][0]-self.beginTime
            t1 = frameInterval[0][1]-self.beginTime
            dt = t1-t0
            locIdx = self.locationIdx[i]
            xOff = self.tcsDict['raOffset'][locIdx]
            yOff = self.tcsDict['decOffset'][locIdx]
            meanTF = self.getTimeBySeq(i)
            fI = {"iframe":i,
                  "begin":t0,
                  "expTime":dt,
                  "loc":locIdx,
                  "offsRA":xOff,
                  "offsDec":yOff,
                  "meanTime":meanTF}
                 #"obsFile":ofs.fileNames[i].obs(),
                 #"ob":ofs.obsFiles[i]}
            if printLines:
                print fI    
            frameInfo.append(fI)
        
        self.frameDict = frameInfo
        return frameInfo
    

    def getSpectralCubeByFrame(self,iFrame,weighted=False, fluxWeighted=False,
                                   wvlStart=None,wvlStop=None,
                                   wvlBinWidth=None,energyBinWidth=None,
                                   wvlBinEdges=None,timeSpacingCut=None):
        """
        return the spectral cube for this frame

        call ObsFile.getSpectralCube for each ObsFile in this frame.
        The dictionary returned copies 'wvlBinEdges' from the first ObsFile,
        and sums the 'cube' and 'effIntTime' from all ObsFiles.

        I left the print statements in to report progress, becuase this is very slow.

        """
        retval = None
        thisInterval = self.frameIntervals[iFrame]
        for i,ofInterval in enumerate(self.obsIntervals):
            overlap = thisInterval & ofInterval
            if len(overlap) > 0:
                tBeg = overlap[0][0]
                tEnd = overlap[0][1]
                integrationTime = tEnd-tBeg
                firstSec = tBeg - self.obsFiles[i].getFromHeader('unixtime')
                obs = self.obsFiles[i]
                obs.setWvlCutoffs(wvlLowerLimit=wvlStart, wvlUpperLimit=wvlStop)
                print time.strftime("%c"),"now call getSpectralCube:  firstSec=",firstSec," integrationTime=",integrationTime, "weighted=",weighted
                spectralCube = obs.getSpectralCube(firstSec=firstSec,
                                                   integrationTime=integrationTime,
                                                   weighted=weighted,
                                                   fluxWeighted=fluxWeighted,
                                                   wvlStart = wvlStart,
                                                   wvlStop = wvlStop,
                                                   wvlBinWidth=wvlBinWidth,
                                                   energyBinWidth=energyBinWidth,
                                                   wvlBinEdges=wvlBinEdges,
                                                   timeSpacingCut=timeSpacingCut
                                                   )
                cube = spectralCube['cube']
                wbe = spectralCube['wvlBinEdges']
                eit = spectralCube['effIntTime']
                if retval is None:
                    retval = {'cube':cube, 'wvlBinEdges':wbe, 'effIntTime':eit}
                else:
                    retval['cube'] += cube
                    retval['effIntTime'] += eit                
        return retval

    def loadSpectralCubes(self,weighted=True, fluxWeighted=False,
                          wvlStart=None,wvlStop=None,
                          wvlBinWidth=None,energyBinWidth=None,
                          wvlBinEdges=None,timeSpacingCut=None):
        """
        calls getSpectralCubeByFrame on each iFrame, storing the
        resuls in the list self.cubes

        use a pickle file named name.pkl as a buffer.  If that file
        exists, load the cubes from there, and save the cubes there
        after loading.  WARNING -- the settings are not stored, so
        they are ignored when loading from the pickle file.

        """
        cpfn = self.name+"-cubes.pkl"
        if os.path.isfile(cpfn):
            print "loadSpectralCubes:  load from ",cpfn
            self.cubes = pickle.load(open(cpfn,'rb'))
        else:
            self.cubes = []
            for iFrame in range(len(self.frameIntervals)):
                print "now load spectral cube for iFrame=",iFrame
                self.cubes.append(self.getSpectralCubeByFrame(iFrame,
                                                              weighted,
                                                              fluxWeighted,
                                                              wvlStart,
                                                              wvlStop,
                                                              wvlBinWidth,
                                                              energyBinWidth,
                                                              wvlBinEdges,
                                                              timeSpacingCut))
            pickle.dump(self.cubes,open(cpfn,'wb'))
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
        
    def loadImageStack(self,fileName, wvlStart=None,wvlStop=None,
                           weighted=False, fluxWeighted=False, 
                           getRawCount=False, scaleByEffInt=False,
                           deadTime=100.e-6):
        #If the file exists, read it out
        if os.path.isfile(fileName):
            print 'Loading image stack from ',fileName
            stackFile = tables.openFile(fullFilename, mode='r')
            images = stackFile.getNode('/',imagesGroupName)._f_getChild(imagesTableName).read()
            images = np.rollaxis(images,2,0)
            try:
                pixIntTimes = stackFile.getNode('/',imagesGroupName)._f_getChild(pixIntTimeTableName).read()
                pixIntTimes = np.rollaxis(pixIntTimes,2,0)
            except tables.exceptions.NoSuchNodeError:
                pixIntTimes = None
            startTimes = stackFile.getNode('/',imagesGroupName)._f_getChild(timeTableName).read()
            try:
                endTimes = stackFile.getNode('/',imagesGroupName)._f_getChild(endTimeTableName).read()
            except tables.exceptions.NoSuchNodeError:
                endTimes = None
            try:
                intTimes = stackFile.getNode('/',imagesGroupName)._f_getChild(intTimeTableName).read()
            except tables.exceptions.NoSuchNodeError:
                intTimes=None
            stackFile.close()
            return {'images':images,'pixIntTimes':pixIntTimes,'startTimes':np.asarray(startTimes).flatten(),'endTimes':np.asarray(endTimes).flatten(),'intTimes':np.asarray(intTimes).flatten()}
        #if the file doesn't exists, make it
        else:
            images=[]
            pixIntTimes=[]
            startTimes=[]
            endTimes=[]
            intTimes=[]
            for iFrame in range(len(self.frameIntervals)):
                im_dict = getPixelCountImageByFrame(self,iFrame,wvlStart,wvlStop,
                                                    weighted, fluxWeighted, 
                                                    getRawCount, scaleByEffInt,
                                                    deadTime)
                images.append(im_dict['images'])
                pixIntTimes.append(im_dict['pixIntTimes'])
                startTimes.append(im_dict['startTimes'])
                endTimes.append(im_dict['endTimes'])
                intTimes.append(im_dict['intTimes'])
                
            
            
            return {'images':images,'pixIntTimes':pixIntTimes,'startTimes':startTimes,'endTimes':endTimes,'intTimes':intTimes}

    def getPixelCountImageByFrame(self,iFrame,wvlStart=None,wvlStop=None,
                           weighted=False, fluxWeighted=False, 
                           getRawCount=False, scaleByEffInt=False,
                           deadTime=100.e-6):
        '''
        This gets the i'th image
        
        Inputs:
            iFrame - which frame you want
            wvlStart - starting wavelength range
            wvlStop - ending wavelength range
            weighted, fluxWeighted, getRawCount, scaleByEffInt - options for obsFile.getPixelCountImage()
            deadTime - for deadtime correcting image
            
        Returns:
            Dictionary with the following keys:
            'image' - fully calibrated, corrected image. 
                      scaled to the total integration time
                      deadTime corrected
            'pixIntTime' - actual integration time for each pixel in image
            'intTime' - length of exposure
            'startTime' - beginning of image (unix time)
            'endTIme' - end of image. Might be different from startTime+intTime if there's a break in the middle while switching to a new obsFile
            
        '''
                           
        scaleByEffInt=False #Do this manually so we can correct deadTime first                           
        retval = None
        for obsInfo in self.frameObsInfos[iFrame]:
            if len(overlap) > 0:
                print obsInfo
                obsInfo['obs'].setWvlCutoffs(wvlLowerLimit=wvlStart, wvlUpperLimit=wvlStop)
                im_dict = obsInfo['obs'].getPixelCountImage(obsInfo["firstSec"],
                                                        obsInfo["integrationTime"],
                                                        weighted,
                                                        fluxWeighted,
                                                        getRawCount,
                                                        scaleByEffInt)
                im = im_dict['image']
                #Correct for dead time
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore",'invalid value encountered in divide',RuntimeWarning)
                    #warnings.simplefilter("ignore",RuntimeWarning)
                    #warnings.simplefilter("ignore",FutureWarning)
                    w_deadTime = 1.0-im_dict['rawCounts']*deadTime/im_dict['effIntTimes']
                im = im/w_deadTime
                #Correct for exposure time
                im = im*obsInfo["integrationTime"]/im_dict['effIntTimes']
                #Remove any funny values
                im[np.invert(np.isfinite(im))]=0.
                
                if retval is None:
                    retval = {'image':im, 'pixIntTime':im_dict['effIntTimes'], 'intTime':obsInfo["integrationTime"],'startTime':self.frameIntervals[iFrame][0],'endTime':self.frameIntervals[iFrame][1]}
                else:
                    retval['image'] += im
                    retval['pixIntTime'] += im_dict['effIntTimes']
                    retval['intTime'] += obsInfo["integrationTime"]
                
        return retval

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


#added by Neil vvvvvv 2/3/2015

class MouseMonitor():
    def __init__(self):
        pass

    def on_click(self,event):
        if event.inaxes is self.ax1:
            self.xyguess1 = [event.xdata,event.ydata]
            print 'Clicked: ',self.xyguess1
        
        elif event.inaxes is self.ax2:
            self.xyguess2 = [event.xdata,event.ydata]
            print 'Clicked: ',self.xyguess2
    
    def on_scroll_cbar(self,event):
        if event.inaxes is self.fig1.cbar.ax:
            increment=0.05
            currentClim = self.fig1.cbar.mappable.get_clim()
            currentRange = currentClim[1]-currentClim[0]
            if event.button == 'up':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]+increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]+increment*currentRange)
            if event.button == 'down':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]-increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]-increment*currentRange)
            self.fig1.cbar.mappable.set_clim(newClim)
            self.fig1.canvas.draw()
           
        elif event.inaxes is self.fig2.cbar.ax:
            increment=0.05
            currentClim = self.fig2.cbar.mappable.get_clim()
            currentRange = currentClim[1]-currentClim[0]
            if event.button == 'up':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]+increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]+increment*currentRange)
            if event.button == 'down':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]-increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]-increment*currentRange)
            self.fig2.cbar.mappable.set_clim(newClim)
            self.fig2.canvas.draw()
    
    def connect(self):
        self.cid1 = self.fig1.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig1.cbar = self.fig1.colorbar(self.handleMatshow1)
        cid1 = self.fig1.canvas.mpl_connect('scroll_event', self.on_scroll_cbar)

        self.cid2 = self.fig2.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig2.cbar = self.fig2.colorbar(self.handleMatshow2)
        cid2 = self.fig2.canvas.mpl_connect('scroll_event', self.on_scroll_cbar)


def getUserObjectGuess(images, norm = None):
    '''
    This is designed to allow the user to look at two frames, determine whether they contain matching stars,
    and if so click on their position withing the image.
    '''
    flagList = np.array([1,1])
    xyguess1 = [0,0]
    xyguess2 = [0,0]

    image1 = images[0]
    image2 = images[1]

    map = MouseMonitor()
    map.fig1 = plt.figure(1)
    map.ax1 = map.fig1.add_subplot(111)
    map.ax1.set_title('Star Position Guess')
    map.handleMatshow1 = map.ax1.matshow(image1,cmap = mpl.cm.gnuplot2, origin = 'lower', norm=norm)
    
    map.fig2 = plt.figure(2)
    map.ax2 = map.fig2.add_subplot(111)
    map.ax2.set_title('Star Position Guess')
    map.handleMatshow2 = map.ax2.matshow(image2,cmap = mpl.cm.gnuplot2, origin = 'lower', norm=norm)

    map.connect()
    
    plt.show() 
    
    try:
        xyguess1 = map.xyguess1
        print 'Guess1 = ' + str(xyguess1)
        flagList[0]=0
    except AttributeError:
        pass

    try:     
        xyguess2 = map.xyguess2
        print 'Guess2 = ' + str(xyguess2)
        flagList[1]=0
    except AttributeError:
        pass
    
    xyguesses = np.array([xyguess1,xyguess2])    

    return xyguesses, flagList

def ObjectFinder(frames, data, RA = None, Dec = None, radiusOfSearch=10, usePsfFit=False):
    '''
    This allows the user to determine the exact positioning of a star withing a given frame, given that
    the star exists in two frames. two images pop up, if they contain the same star, click on it in both images,
    and the exact pixel positioning will be determine via centroiding. Some parts are fudged for now.

    Inputs:
        frames - This is an array of images corresponding to the various frames used for the mosaic. note that 
                 the way this is written as of now requires that the first  and last frame include the same star. 
                 I think this is typical so it should not be a problem.

        data - this is an array of dictionaries corresponding to each frame. This is generated from
               getFrameDict() which I have included in the ObsFileSeq class. This is required.

        RA/Dec - these are not used as of now.
    
    Outputs:
        
        all the parameters that are needed for setRm()

    Ben also suggested that a cross-correlating technique may be useful for matching frames with no star in them.
    What do you guys think? I need to look into it more but I believe that this code can be expanded to also do
    cross-correlation. - Neil
    '''

    offsRA = []
    offsDec = []
    iframe = [] 
    meanTime = []  
    
    for i in range(len(data)):
        offsRA.append(data[i]["offsRA"])
        offsDec.append(data[i]["offsDec"])
        iframe.append(data[i]["iframe"])
        meanTime.append(data[i]["meanTime"])

    offsRA = np.array(offsRA)
    offsDec = np.array(offsDec)
    iframe = np.array(iframe)
    meanTime = np.array(meanTime)
    
    dpp = []
    ang = []
    
    for i in range(1, len(frames)):
        images = np.array([frames[i-1], frames[i]])    
        oRA = np.array([offsRA[i-1], offsRA[i]])
        oDec = np.array([offsDec[i-1], offsDec[i]])
        print 'Looking for matching Stars...'

        print 'frame: ', iframe[i-1], 'Offset RA: ', offsRA[i-1], 'Offset Dec: ', offsDec[i-1]
        print 'frame: ', iframe[i], 'Offset RA: ', offsRA[i], 'Offset Dec: ', offsDec[i]         
        
        xyguesses, flagList = getUserObjectGuess(images)

        if flagList[1]==0 and flagList[0]==0:
            print 'match found! - determining centroid positions'
                
            xycenter1, flag1 = cc.centroidImage(images[1], xyguesses[1], radiusOfSearch=radiusOfSearch, doDS9=False, usePsfFit=usePsfFit)
            xycenter0, flag0 = cc.centroidImage(images[0], xyguesses[0], radiusOfSearch=radiusOfSearch, doDS9=False, usePsfFit=usePsfFit)
            print 'Success! Matching stars at: ', xycenter0, xycenter1                
                
            rc1 = np.array(xycenter1)
            rc0 = np.array(xycenter0)
            dCol = rc1[0] - rc0[0]
            dRow = rc1[1] - rc0[1]                
            dPix = math.sqrt(((rc1-rc0)**2).sum())

            #center ra,dec of fram calculated from offsets
            dRA = oRA[1]-oRA[0] #arcseconds
            dDec = oDec[1]-oDec[0]
            dDeg = math.sqrt(dRA**2+dDec**2)/3600 #degrees
                
            degPerPix = dDeg/dPix #plate scale
    
            #rotation
            thetaPix = math.atan2(dCol,dRow) #angle from verticle
            thetaSky = math.atan2(dRA, dDec) #angle from north
            theta = thetaPix-thetaSky        #degrees
                
            dpp.append(degPerPix)
            ang.append(theta)
        
        elif flagList[1]==1 or flagList[0]==1:
            print 'no star found' 

    dpp = np.array(dpp)
    #print dpp
    degPerPix = np.mean(dpp)
    #print degPerPix
    ang = np.array(ang)
    #print ang
    theta = np.mean(ang)  
    #print theta  
     
    ## Pick two frames where the ra,dec offset is zero,
    # usually the beginning and ending frames
    
    print 'Matching stars from the first and last frames'    
    
    images = [frames[0], frames[-1]]

    print 'frame: ', iframe[0], 'Offset RA: ', offsRA[0], 'Offset Dec: ', offsDec[0]
    print 'frame: ', iframe[-1], 'Offset RA: ', offsRA[-1], 'Offset Dec: ', offsDec[-1]  
    
    xyguesses, flagList = getUserObjectGuess(images)
    
    xycenter1, flag1 = cc.centroidImage(images[1], xyguesses[1], radiusOfSearch=radiusOfSearch, doDS9=False, usePsfFit=usePsfFit)
    xycenter0, flag0 = cc.centroidImage(images[0], xyguesses[0], radiusOfSearch=radiusOfSearch, doDS9=False, usePsfFit=usePsfFit)
    
    print 'Success! Matching stars at: ', xycenter0, xycenter1        

    rcA = np.array(xycenter1)
    rcB = np.array(xycenter0)
    sct = math.cos(theta)*degPerPix    
    sst = math.sin(theta)*degPerPix
    # This rotation matrix converts from row,col to ra,dec in degrees
    rm = np.array([[sct,-sst],[sst,sct]])
    rdA = rm.dot(rcA)
    rdB = rm.dot(rcB)
    deltaRa = rdB[0]-rdA[0]
    deltaTime = meanTime[-1] - meanTime[0]
    raArcsecPerSec = 3600*deltaRa/deltaTime
    
    return degPerPix, theta, raArcsecPerSec
       

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
    print "now get time of first frame"
    for i in range(66):
        print "i="," time=",ofs.getTimeBySeq(i)
    #apci = ofs.getAllPixelCountImages(getRawCount=True)
    del ofs
