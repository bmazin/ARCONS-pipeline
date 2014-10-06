"""
Author: Alex Walter
Date: Oct 3, 2014

This is a wrapper class for appyling time masks. The main advantage is the 'mask' keyword
parameter which allows you to mask out time intervals which are masked for a certain reason.
"""





import numpy as np
from interval import interval


class hotPixelMasker:

    def __init__(self, intervals, reasons, reasonEnum, nRow, nCol, obsFileName, ticksPerSec, expTime, startTime, endTime, mask=None):
        """
        Use readHotPixels() in hotPixels.py to create a hotPixelMasker object. 
        
        keywords:
            mask - A list of integers corresponding to a reason enum type in headers/TimeMask.py.
                   If the reason is indicated it will use it 
        """
        
        self.intervals = intervals
        self.reasons=reasons
        self.reasonEnum = reasonEnum
        self.nRow = nRow
        self.nCol = nCol
        self.obsFileName = obsFileName
        self.ticksPerSec = ticksPerSec
        self.expTime = expTime
        self.startTime=startTime
        self.endTime=endTime
        
        if mask==None or mask==-1 or mask=='all':
            self.mask = reasonEnum.__dict__['_names'].values()
        else: self.mask=mask
        
    def set_mask(self,reasons=[]):
        if not len(reasons==0):
            self.mask = [self.reasonEnum[reason] for reason in reasons]
        
    def get_intervals(self,row,col,reasons=[]):
        mask = self.mask if len(reasons)==0 else [self.reasonEnum[reason] for reason in reasons]
        intervalList = [self.intervals[row,col][i] for i in np.where(np.in1d(self.reasons[row,col],mask))[0]]
        return interval.union(intervalList)


    def getEffIntTimeImage(self,firstSec=0,integrationTime=-1,reasons=[]):
    '''    
    Get the total effective exposure time for each pixel after subtracting 
    any intervals where a pixel was masked as hot or bad.
    
    INPUTS:
        firstSec - Start time (sec) to start calculations from, starting from
                  the beginning of the exposure to which timeMask refers.
        integrationTime - Length of integration time (sec) from firstSec to include
                  in the calculation. 
        reasons - Reasons to include in the bad pixel mask. By default, uses whatever
                  is in self.mask

    RETURNS:
        A 2D array representing the total effective exposure time
        in seconds for each pixel in the detector array.
    '''
    
        lastSec = (firstSec + integrationTime) if (integrationTime>0 and (firstSec + integrationTime) < self.expTime) else self.expTime
        integrationInterval = interval([firstSec, lastSec])
        effectiveIntTimes = np.zeros((self.nRow,self.nCol),dtype=float)
        effectiveIntTimes.fill(np.nan)
        
        for iRow in np.arange(self.nRow):
            for iCol in np.arange(self.nCol):
                #Get the unioned (possibly multi-component) bad interval for this pixel.
                #(As in ObsFile.getPixelBadTimes)
                allBadIntervals = self.get_intervals(iRow, iCol,reasons)
                #Get intersection of integration time interval and the bad time intervals.
                maskedIntervals = allBadIntervals & integrationInterval
                effectiveIntTimes[iRow,iCol] = (lastSec-firstSec) - utils.intervalSize(maskedIntervals)
        
        return effectiveIntTimes
    
    def getHotPixels(self,firstSec=0,integrationTime=-1,reasons=[]):
    '''
    Return a boolean array indicating which pixels went bad at any point
    during the specified integration time.
    
    INPUTS:
        firstSec - Start time (sec) to start calculations from, starting from
                    the beginning of the exposure to which timeMask refers.
        integrationTime - Length of integration time (sec) from firstSec to include
                    in the calculation. 
    RETURNS:
        A 2D Boolean array matching the size/shape of the detector image. True indicates
        a pixel that went bad between firstSec and firstSec+integrationTime, and False 
        indicates that the pixel was okay during that time.
    '''
    
        effectiveIntTimes=self.getEffIntTimeImage(firstSec=firstSec,integrationTime=integrationTime,reasons=reasons)
        lastSec = (firstSec + integrationTime) if (integrationTime>0 and (firstSec + integrationTime) < self.expTime) else self.expTime
        return effectiveIntTimes<(lastSec-firstSec)
    




