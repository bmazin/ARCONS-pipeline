from fluxcal.fluxCal import FluxCal
from util.ObsFile import ObsFile
from util.FileName import FileName
import hotpix.hotPixels as hp
#from util import popup
from util.rebin import rebin
from util import MKIDStd
from util import utils
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

"""
Created 2/7/2013 by Seth Meeker
Routine for anaylizng Landolt standard through V and R band filters to check total throughput

"""
def getTimeMaskFileName(obsFileName):
    scratchDir = os.getenv('MKID_PROC_PATH')
    hotPixDir = os.path.join(scratchDir,'timeMasks')
    fileName = obsFileName.split('/')[-1]
    fileNameBase = fileName.split('_')[-1]
    newName = 'timeMask_'+fileNameBase
    fullHotPixFileName = os.path.join(hotPixDir,newName)
    return fullHotPixFileName

def loadSpectra(fileObj,nCol,nRow):
    Spectra = [[[] for i in xrange(nCol)] for j in xrange(nRow)]
    #print fileObj
    for iRow in xrange(nRow):
        for iCol in xrange(nCol):
            #print iRow,iCol,
            count = fileObj.getPixelCount(iRow,iCol)
            #print count
            Spectra[iRow][iCol],wvlBinEdges = fileObj.getPixelSpectrum(iRow,iCol,weighted=True,fluxWeighted=True,firstSec=0,integrationTime=-1)
    Spectra = np.array(Spectra)
    return wvlBinEdges,Spectra

def loadStd(objectName,wvlBinEdges):
    #import the known spectrum of the calibrator and rebin to the histogram parameters given
    #must be imported into array with dtype float so division later does not have error
    std = MKIDStd.MKIDStd()
    a = std.load(objectName)
    a = std.normalizeFlux(a)
    x = a[:,0]
    y = a[:,1]
    realSpectra = y
    realSpectraWvl = x
    newa = rebin(x,y,wvlBinEdges)
    binnedSpectra = newa[:,1]
    return binnedSpectra

def calculateMedian(spectra, nCol, nRow, nWvlBins):
    """given array of spectra (one for each pixl in nCol by nRow grid) returns median spectrum"""
    spectra2d = np.reshape(spectra,[nRow*nCol,nWvlBins])
    wvlMedian = np.empty(nWvlBins,dtype=float)
    for iWvl in xrange(nWvlBins):
        spectrum = spectra2d[:,iWvl]
        goodSpectrum = spectrum[spectrum != 0]#dead pixels need to be taken out before calculating medians
        wvlMedian[iWvl] = np.median(goodSpectrum)
    return wvlMedian

def main():
    """
    params = []
    paramfile = sys.argv[1]

    f = open(paramfile,'r')
    for line in f:
        params.append(line)
    f.close()
    datadir = params[0].split('=')[1].strip()
    flatdir = params[1].split('=')[1].strip()
    fluxdir = params[2].split('=')[1].strip()
    wvldir = params[3].split('=')[1].strip()
    obsfile = params[4].split('=')[1].strip()
    skyfile = params[5].split('=')[1].strip()
    flatfile = params[6].split('=')[1].strip()
    fluxfile = params[7].split('=')[1].strip()
    wvlfile = params[8].split('=')[1].strip()
    objectName = params[9].split('=')[1].strip()
    fluxCalObject = params[10].split('=')[1].strip()

    obsFileName = os.path.join(datadir, obsfile)
    skyFileName = os.path.join(datadir, skyfile)
    wvlCalFileName = os.path.join(wvldir, wvlfile)
    flatCalFileName = os.path.join(flatdir, flatfile)
    fluxCalFileName = os.path.join(fluxdir, fluxfile)
    """

    if len(sys.argv) >2:
        fileNum = str(sys.argv[2])
    else:
        fileNum = '0'

    #science object parameter file
    params = []
    paramfile = sys.argv[1]
    f = open(paramfile,'r')
    for line in f:
        params.append(line)
    f.close()

    datadir = params[0].split('=')[1].strip()
    flatdir = params[1].split('=')[1].strip()
    wvldir = params[2].split('=')[1].strip()
    obsfile = params[3].split('=')[1].strip()
    skyfile = params[4].split('=')[1].strip()
    flatfile = params[5].split('=')[1].strip()
    wvlfile = params[6].split('=')[1].strip()
    objectName = params[9].split('=')[1].strip()

    wvldir = "/Scratch/waveCalSolnFiles/oldbox_numbers/20121206"

    if len(params)>10:
        xpix = int(params[10].split('=')[1].strip())
        ypix = int(params[11].split('=')[1].strip())
        apertureRadius = int(params[12].split('=')[1].strip())
        startTime = int(params[13].split('=')[1].strip())
        intTime =int(params[14].split('=')[1].strip())

    obsFileName = os.path.join(datadir, obsfile)
    skyFileName = os.path.join(datadir, skyfile)
    wvlCalFileName = os.path.join(wvldir, wvlfile)
    flatCalFileName = os.path.join(flatdir, flatfile)

    obs = ObsFile(obsFileName)
    obs.loadWvlCalFile(wvlCalFileName)
    obs.loadFlatCalFile(flatCalFileName)
    print "analyzing file %s"%(obsFileName)
    print "loaded data file and calibrations\n---------------------\n"

    nRow = obs.nRow
    nCol = obs.nCol
    obsTime = obs.getFromHeader("exptime")
    #wvlBinEdges,obsSpectra = loadSpectra(obs,nCol,nRow)
    #nWvlBins=len(wvlBinEdges)-1

    #print np.shape(obsSpectra)
    #print nRow
    #print nCol
    #print nWvlBins

    #load/generate hot pixel mask file
    HotPixFile = getTimeMaskFileName(obsFileName)
    if not os.path.exists(HotPixFile):
        hp.findHotPixels(obsFileName,HotPixFile)
        print "Flux file pixel mask saved to %s"%(HotPixFile)
    obs.loadHotPixCalFile(HotPixFile)
    print "Hot pixel mask loaded %s"%(HotPixFile)

    #######
    #EVERYTHING BEFORE HERE IS STANDARD FILE/CALFILE LOADING

    
    startWvl = 3000
    #stopWvl = 7000 #for V-band
    stopWvl = 9000 #for R-band

    print "Making spectral cube"
    #for pg0220 first sec should be 80 since object is moving around before this
    #for pg0220A first sec should be 70, integration time is 140
    #for landolt 9542 first sec should be 20, int time is -1
    cubeDict = obs.getSpectralCube(firstSec=startTime, integrationTime=intTime, wvlStart = startWvl, wvlStop = stopWvl, wvlBinEdges = [startWvl,stopWvl], weighted=False)

    cube= np.array(cubeDict['cube'], dtype=np.double)
    wvlBinEdges= cubeDict['wvlBinEdges']
    effIntTime= cubeDict['effIntTime']
    print "median effective integration time = ", np.median(effIntTime)

    nWvlBins=len(wvlBinEdges)-1
    print "cube shape ", np.shape(cube)
    print "effIntTime shape ", np.shape(effIntTime)

    #add third dimension to effIntTime for  broadcasting
    effIntTime = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
    cube /= effIntTime #put cube into counts/s
    
    #calculate midpoints of wvl bins for plotting
    wvls = np.empty((nWvlBins),dtype=float)
    for n in xrange(nWvlBins):
        binsize=wvlBinEdges[n+1]-wvlBinEdges[n]
        wvls[n] = (wvlBinEdges[n]+(binsize/2.0))

    print "wvls ",wvls
    #reshape cube for makeMovie
    movieCube = np.zeros((nWvlBins,np.shape(cube)[0],np.shape(cube)[1]),dtype=float)
    for i in xrange(nWvlBins):
        movieCube[i,:,:] = cube[:,:,i]
        #show individual frames as they are made to debug
        #plt.matshow(movieCube[i],vmin = 0, vmax = 100)
        #plt.show()
    print "movieCube shape ", np.shape(movieCube)
    print "wvls shape ", np.shape(wvls)

    #print cube
    #print "--------------------------"
    #print movieCube

    print "adding frames with wvl below ", stopWvl
    finalframe = np.zeros((1,np.shape(movieCube)[1],np.shape(movieCube)[2]))
    for f in xrange(len(wvls[wvls<stopWvl])):
        print wvls[f]
        finalframe[0]+=movieCube[f]
        plt.matshow(movieCube[f],vmin=0,vmax = 40)
        plt.show()

    movieCube = finalframe

    np.savez('%s_%s.npz'%(objectName,fileNum),stack=movieCube,wvls=wvls)
    print "Saved frame to .npz file"
    
    plt.matshow(movieCube[0],vmin=0,vmax = 40)
    plt.show()

if __name__ == '__main__':
    main()














