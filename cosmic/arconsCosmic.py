#
# Look at a .h5 file from the Lick 2012 run
# Also, make a histogram of all the peak fits, which cna then be converted to an energy

# Use the /beammap/beamimage information to get the list of roach board/pixel/time
#
# Set the environment variable MKID_RAW_PATH to point to the data location
#
# This is the format of a data packet:
#  8 bits for channel
# 12 bits for parabola fit peak
# 12 bits for data peak
# 12 bits for baseline
# 20 bits for timestamp

import sys, os
import tables
import pyfits
import numpy as np
import matplotlib.pyplot as plt

def countPhotons(hdf5FileName):
    # make the full file name by joining the input name to the MKID_RAW_PATH (or .)
    dataDir = os.getenv('MKID_RAW_PATH','.')
    hdf5FullFileName = os.path.join(dataDir,hdf5FileName)
    fid = tables.openFile(hdf5FullFileName, mode='r')
    header = fid.root.header.header
    titles = header.colnames
    info = header[0]
    exptime = info[titles.index('exptime')]

    timeMask = int(20*'1',2)   #bitmask of 20 ones
    timeMin = timeMask
    timeMax = 0
    nPhoton = 0

    print "before loadAllSecs"
    nRow,nCol,allSecs,beamImage = loadAllSecs(fid)
    print "after  loadAllSecs"
    nPhoton = 0
    for iRow in range(nRow):
        for iCol in range(nCol):
            for sec in allSecs[iRow,iCol]:
                nPhoton += len(sec)
    fid.close()
    return nPhoton

def getHgByRoach(hdf5FileName):
    hgs = {}
    # make the full file name by joining the input name to the MKID_RAW_PATH (or .)
    dataDir = os.getenv('MKID_RAW_PATH','.')
    hdf5FullFileName = os.path.join(dataDir,hdf5FileName)
    fid = tables.openFile(hdf5FullFileName, mode='r')
    header = fid.root.header.header
    titles = header.colnames
    info = header[0]
    exptime = info[titles.index('exptime')]

    timeMask = int(20*'1',2)   #bitmask of 20 ones
    timeMin = timeMask
    timeMax = 0
    nPhoton = 0

    print "before loadAllSecs"
    nRow,nCol,allSecs,beamImage = loadAllSecs(fid)
    print "after  loadAllSecs"
    nPhoton = 0
    nBins = 10
    for iSec in range(exptime):
        print "iSec=%4d / %4d" % (iSec,exptime)
        hgsThisSec = {}
        for iRow in range(nRow):
            for iCol in range(nCol):
                sec = allSecs[iRow,iCol][iSec]
                if len(sec) > 0:
                    times = sec & timeMask
                    hg,edges = np.histogram(times,bins=nBins,range=(0,10000))
                    roachName = beamImage[iRow][iCol].split("/")[1]
                    if not hgsThisSec.has_key(roachName):
                        hgsThisSec[roachName] = np.zeros(nBins,dtype=np.int64)
                    hgsThisSec[roachName] += hg
                    nPhoton += len(sec)
        for roachName in hgsThisSec.keys():
            if not hgs.has_key(roachName):
                hgs[roachName] = []
            hgs[roachName] += list(hgsThisSec[roachName])
    fid.close()
    times = np.arange(0,1,step=1.0/nBins)
    return nPhoton,hgs,times

def plotHgByRoach(hdf5FileName,hgs,times):
    pfn = os.path.splitext(hdf5FileName)[0]+"-hgByRoach.png"
    plt.clf()
    keys = hgs.keys()
    keys.sort()
    for roachName in keys:
        print "roachName=",roachName
        plt.plot(times,hgs[roachName],label=roachName)
    plt.title(hdf5FileName)
    dt = times[1]-times[0]
    ylabel = "counts per %.1f sec" % dt
    plt.ylabel(ylabel)
    plt.xlabel("time (seconds)")
    plt.xlim(180,260)
    plt.legend(loc='best', fancybox=True)
    plt.legend().get_frame().set_alpha(0.5)
    #plt.savefig(pfn)
    

def getMinMaxTime(hdf5FileName):
    # make the full file name by joining the input name to the MKID_RAW_PATH (or .)
    dataDir = os.getenv('MKID_RAW_PATH','.')
    hdf5FullFileName = os.path.join(dataDir,hdf5FileName)
    fid = tables.openFile(hdf5FullFileName, mode='r')
    header = fid.root.header.header
    titles = header.colnames
    info = header[0]
    exptime = info[titles.index('exptime')]

    timeMask = int(20*'1',2)   #bitmask of 20 ones
    timeMin = timeMask
    timeMax = 0
    nPhoton = 0

    print "before loadAllSecs"
    nRow,nCol,allSecs,beamImage = loadAllSecs(fid)
    print "after  loadAllSecs"
    bins = 100
    hgSum = np.zeros(bins)
    for iSecond in np.arange(exptime):
        for iRow in np.arange(nRow):
            for iCol in np.arange(nCol):
                secs = allSecs[iRow,iCol]
                print "iRow=%02d iCol=%02d pixel=%20s  nSecs=%5d %10d %10d" % (iRow, iCol, pixel, len(secs), timeMin, timeMax)
                sec = secs[iSecond]
                if len(sec) > 0:
                    nPhoton += len(sec)
                    times = sec & timeMask
                    hg,edges = np.histogram(times,bins=bins,range=(0,1000000))
                    hgSum += hg
                    timeMin = min(times.min(),timeMin)
                    timeMax = max(times.max(),timeMax)
    fid.close()
    return timeMin,timeMax,nPhoton,hgSum

def loadAllSecs(fid):
    beamImage = fid.getNode("/beammap/beamimage")
    nRow = len(beamImage)
    nCol = len(beamImage[0])
    allSecs = dict( ((i,j),None) for i in range(nRow) for j in range(nCol) )
    for iRow in np.arange(nRow):
        for iCol in np.arange(nCol):
            allSecs[iRow,iCol] = fid.getNode(beamImage[iRow][iCol])
    return nRow, nCol, allSecs, beamImage

def makeOrLoadPklFile(hdf5FileName):
    pklFileName = hdf5FileName+".pkl"
    if os.path.exists(pklFileName) :
        print "Read from pklFileName=",pklFileName
        pFile = open(pklFileName,'r')
        nPhoton = pickle.load(pFile)
        pixels = pickle.load(pFile)
        print "  ..parabolaFitPeaks"
        parabolaFitPeaks = pickle.load(pFile)
        print "  ..dataPeaks"
        dataPeaks = pickle.load(pFile)
        print "  ..baselines"
        baselines = pickle.load(pFile)
        print "  ..delaTimes"
        deltaTimes = pickle.load(pFile)
        return nPhoton, pixels, parabolaFitPeaks, dataPeaks, baselines, deltaTimes

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
    nRow = shape[0]
    nCol = shape[1]
    parabolaFitPeaks = np.zeros((nRow,nCol,2**12), dtype=np.uint32)
    dataPeaks = np.zeros((nRow,nCol,2**12), dtype=np.uint32)
    baselines = np.zeros((nRow,nCol,2**12), dtype=np.uint32)
    deltaTimeMax=2**10
    deltaTimes = np.zeros((nRow,nCol,deltaTimeMax), dtype=np.uint32)
    nPhoton = 0
    for iRow in np.arange(len(beamImage)):
        for iCol in np.arange(len(beamImage[0])):
            if (iRow>=0 and iCol>=0):
                print "gather photons from iRow=",iRow," iCol=", iCol
                pixel = beamImage[iRow][iCol]                
                secs = fid.getNode("/"+pixel)
                channels, sum, hgParabolaFitPeak, hgDataPeak, hgBaseline, hgDeltaTime  = iterateOverPixel(secs)
                print "nPhoton=",nPhoton,"   channels=",channels
                nPhoton += sum
                pixels[iRow][iCol] = sum
                parabolaFitPeaks[iRow,iCol,:] = hgParabolaFitPeak
                dataPeaks[iRow,iCol,:] = hgDataPeak
                baselines[iRow,iCol,:] = hgBaseline
                deltaTimes[iRow,iCol,:] = hgDeltaTime
    print "Write pkl file ",pklFileName
    pFile = open(pklFileName,'wb')
    pickle.dump(nPhoton,pFile)
    pickle.dump(pixels,pFile)
    print "   ..parabolaFitPeaks"
    pickle.dump(parabolaFitPeaks,pFile)
    print "   ..dataPeaks"
    pickle.dump(dataPeaks,pFile)
    print "   ..baseliness"
    pickle.dump(baselines,pFile)
    print "   ..deltaTimes"
    pickle.dump(deltaTimes,pFile)
    print "   ..close"
    pFile.close()

    fid.close()
    return nPhoton, pixels, parabolaFitPeaks, dataPeaks, baselines, deltaTimes


def iterateOverPixel(secs):
    channels = []
    sum = 0
    hgParabolaFitPeak =  np.zeros(2**12, dtype=np.uint32)
    hgDataPeak = np.zeros(2**12, dtype=np.uint32)
    hgBaseline = np.zeros(2**12, dtype=np.uint32)
    hgDeltaTime = np.zeros(2**10, dtype=np.uint32)

    pulseMask = int(12*'1',2)  #bitmask of 12 ones
    timeMask = int(20*'1',2)   #bitmask of 20 ones
    
    for sec in secs:
        previousTime = -9999
        for packet in sec:
            # here is one 64-bit packet.  
            packet = int(packet)
            # unpack the packet with the Lick 2012 format
            channel = packet>>56
            parabolaFitPeak = packet>>44 & pulseMask
            dataPeak = packet>>32 & pulseMask
            baseline = packet>>20 & pulseMask
            time = packet & timeMask

            # Tally the values
            if (channels.count(channel) == 0) :
                channels.append(channel)
            sum += 1
            hgParabolaFitPeak[parabolaFitPeak-baseline] += 1
            hgDataPeak[dataPeak-baseline] += 1
            hgBaseline[baseline] += 1

            if (previousTime != -9999) :
                dt = time - previousTime
                if dt < 0:
                    dt += 2**20
                if (dt >= 0 and dt < hgDeltaTime.size-1):
                    hgDeltaTime[dt] += 1
            previousTime = time
    return channels, sum, hgParabolaFitPeak, hgDataPeak, hgBaseline, hgDeltaTime

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
