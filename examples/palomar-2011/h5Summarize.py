#
# Look at a .h5 file from the Palomar 2011 run
# Also, make a histogram of all the peak fits, which cna then be converted to an energy

# Use the /beammap/beamimage information to get the list of roach board/pixel/time
#
# Set the environment variable MKID_RAW_PATH to point to the data location
#
# Example use:
# 

import sys, os
import tables
import pyfits
import numpy as np
import matplotlib.pyplot as plt
import pickle

def peakfit(y1,y2,y3):
	#returns the y value of the peak of the parabola defined by the given 3 points
	#where 3 points are separated by 1 in x
    y4=y2-0.125*((y3-y1)**2)/(y3+y1-2*y2)
    return y4


def makeOrLoadPklFile(hdf5FileName):
    pklFileName = hdf5FileName+".pkl"
    #if os.path.exists(pklFileName) :
    #    # print("pklFileName="+pklFileName+" already exists")
    #    pFile = open(pklFileName,'r')
    #    nPhoton = pickle.load(pFile)
    #    nDtNegative = pickle.load(pFile)
    #    hgTime = pickle.load(pFile)
    #    hgAt = pickle.load(pFile)
    #    pixels = pickle.load(pFile)
    #    pFile.close()
    #    return nPhoton, nDtNegative, hgTime, hgAt, pixels
    # make the full file name by joining the input name to the MKID_RAW_PATH (or .)
    dataDir = os.getenv('MKID_RAW_PATH','.')
    hdf5FullFileName = os.path.join(dataDir,hdf5FileName)
    print "full file name is ",hdf5FullFileName
    if (not os.path.exists(hdf5FullFileName)):
        print "file does not exist: ",hdf5FullFileName
        sys.exit(1)

    # open the actual file.  This might take a while
    fid = tables.openFile(hdf5FullFileName, mode='r')

    # get the header information
    header = fid.root.header.header
    titles = header.colnames
    info = header[0]
    # get the beam image.  This is a 2d array of roach board/pixel/time locations
    try:
        beamImage = fid.getNode("/beammap/beamimage")
    except:
        print "/beammap/beamimage"
        return 0,0,[],[],[[]]
    # Make a 2d array to hold the image
    shape = beamImage.shape
    pixels = np.zeros(shape,dtype=np.uint32)
    hgTimes = []
    hgAts = []
    # count the total number of photons in the file
    nPhoton = 0
    nDtNegative = 0
    for iRow in np.arange(len(beamImage)):
        for iCol in np.arange(len(beamImage[0])):
            if (iRow>=0 and iCol>=0):
                pixel = beamImage[iRow][iCol]
                secs = fid.getNode(pixel)
                sum, thisNDtNegative, hgAt, hgTime = iterateOverPixel(secs)
                pixels[iRow][iCol] = sum
                hgTimes.append(hgTime)
                hgAts.append(hgAt)
                nDtNegative += thisNDtNegative
                nPhoton += sum

    pFile = open(pklFileName,'wb')
    pickle.dump(nPhoton,pFile)
    pickle.dump(nDtNegative,pFile)
    pickle.dump(hgTimes,pFile)
    pickle.dump(hgAts,pFile)
    pickle.dump(pixels,pFile)
    pFile.close()

    fid.close()
    return nPhoton, nDtNegative, hgTimes, hgAts, pixels

def iterateOverPixel(secs):
    # Make a 1d array to hold the fit peak
    hgFit = np.zeros(1000+2**12, dtype=np.uint32)
    # Make a 1d array to hold the at peak values
    hgAt = np.zeros(1000+2**12, dtype=np.uint32)
    # Make a 1d array to hold the delta time distribution
    hgTime = np.zeros(2**13, dtype=np.uint32)
    pulseMask = int(12*'1',2)  #bitmask of 12 ones
    timeMask = int(20*'1',2)   #bitmask of 20 ones
    sum = 0
    thisNDtNegative = 0
    
    for sec in secs:
        previousTime = -9999
        for packet in sec:
            # here is the 64-bit number.  
            packet = int(packet)
            beforePeak = packet>>44 & pulseMask
            atPeak = packet>>32 & pulseMask
            afterPeak = packet>>20 & pulseMask
            peak = peakfit(beforePeak, atPeak, afterPeak)
            hgFit[peak] += 1
            hgAt[atPeak] += 1
            sum += 1
            time = packet & timeMask
            if (previousTime != -9999) :
                dt = time - previousTime
                if dt < 0:
                    dt += 2**20
                    if (dt < 0):
                        thisNDtNegative += 1
                    elif (dt < hgTime.size-1):
                        hgTime[dt] += 1
            previousTime = time
    return sum, thisNDtNegative, hgTime, hgAt

def makeAllPklFiles():
    dataDir = os.getenv('MKID_RAW_PATH','.')
    for h5File in [f for f in os.listdir(dataDir) if f.endswith(".h5")]:
        makeOrLoadPklFile(h5File)

def makeSummaryTable(stName):

    stFile = open(stName,'w')

    line = "%(fName)25s  %(nPhoton)10s %(nDtNegative)10s %(iMinTime)5s %(nRow)3s %(nCol)3s %(nZeroPix)4s %(target)s" % \
        {'fName': 'h5File', "nPhoton": 'nPho', "nDtNegative":'nDtN', "iMinTime":'minT', \
             'nRow':'nR', 'nCol':'nC', 'nZeroPix':'n0P', 'target':'target'}
    print line
    stFile.write(line+"\n")
    dataDir = os.getenv('MKID_RAW_PATH','.')
    for h5File in [f for f in os.listdir(dataDir) if f.endswith(".h5")]:
        dataDir = os.getenv('MKID_RAW_PATH','.')
        hdf5FullFileName = os.path.join(dataDir,h5File)
        fid = tables.openFile(hdf5FullFileName, mode='r')
        header = fid.root.header.header
        titles = header.colnames
        info = header[0]
        try:
            index = titles.index('target')
            target = info[index]
        except:
            target = 'unknown'

        fid.close()
        nPhoton, nDtNegative, hgTime, hgAt, pixels = makeOrLoadPklFile(h5File)
        iMinTime = -1
        for i in np.arange(len(hgTime)):
            if (hgTime[i] > 0):
                iMinTime = i
                break

        nRow = len(pixels)
        nCol = len(pixels[0])
        nZeroPix = 0
        for row in pixels:
            for v in row:
                if v == 0:
                    nZeroPix += 1

        line = "%(fName)25s  %(nPhoton)10d %(nDtNegative)10d %(iMinTime)5d %(nRow)3d %(nCol)3d %(nZeroPix)4d %(target)s" % \
            {'fName': h5File, "nPhoton": nPhoton, "nDtNegative":nDtNegative, "iMinTime":iMinTime, \
                 'nRow':nRow, 'nCol':nCol, 'nZeroPix':nZeroPix, 'target':target}
        print line
        stFile.write(line+"\n")
    stFile.close()
