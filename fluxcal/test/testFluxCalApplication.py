from fluxcal.fluxCal import FluxCal
from util.ObsFile import ObsFile
from util.FileName import FileName
import hotpix.hotPixels as hp
#from util import popup
from util.rebin import rebin
from util import MKIDStd
import matplotlib.pyplot as plt
import numpy as np
import os, sys
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

"""
Created 2/7/2013 by Seth Meeker
Routine for testing the application of FluxCal files.

"""
def getTimeMaskFileName(obsFileName):
    scratchDir = os.getenv('INTERM_PATH')
    hotPixDir = os.path.join(scratchDir,'timeMasks')
    fileName = obsFileName.split('/')[-1]
    fileNameBase = fileName.split('_')[-1]
    newName = 'timeMask_'+fileNameBase
    fullHotPixFileName = os.path.join(hotPixDir,newName)
    return fullHotPixFileName

def loadSpectra(fileObj,nCol,nRow):
    Spectra = [[[] for i in xrange(nCol)] for j in xrange(nRow)]
    effIntTime = [[[] for i in xrange(nCol)] for j in xrange(nRow)]
    #print fileObj
    for iRow in xrange(nRow):
        for iCol in xrange(nCol):
            #print iRow,iCol,
            count = fileObj.getPixelCount(iRow,iCol)
            #print count
            specDict = fileObj.getPixelSpectrum(iRow,iCol,weighted=True,fluxWeighted=True,firstSec=0,integrationTime=-1)
            Spectra[iRow][iCol],wvlBinEdges,effIntTime[iRow][iCol] = specDict['spectrum'],specDict['wvlBinEdges'],specDict['effIntTime'] 
    Spectra = np.array(Spectra)
    return wvlBinEdges,Spectra,effIntTime

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

    if len(sys.argv) >3:
        filenum = str('_'+sys.argv[3])
    else:
        filenum = '_0'

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

    if len(params)>10:
        xpix = params[10].split('=')[1].strip()
        ypix = params[11].split('=')[1].strip()
        apertureRadius = params[12].split('=')[1].strip()


    #flux cal object parameter file

    params2 = []
    param2file = sys.argv[2]
    f = open(param2file,'r')
    for line in f:
        params2.append(line)
    f.close()

    fluxdir = params2[7].split('=')[1].strip()
    fluxfile = params2[8].split('=')[1].strip()
    fluxCalObject = params2[9].split('=')[1].strip()

    obsFileName = os.path.join(datadir, obsfile)
    skyFileName = os.path.join(datadir, skyfile)
    wvlCalFileName = os.path.join(wvldir, wvlfile)
    flatCalFileName = os.path.join(flatdir, flatfile)
    fluxCalFileName = os.path.join(fluxdir, fluxfile)

    print "obsfile = ",obsFileName
    print "skyfile = ",skyFileName
    print "wvlcal = ", wvlCalFileName
    print "flatcal = ", flatCalFileName
    print "fluxcal = ", fluxCalFileName
    print "object = ", objectName
    print "flux cal object = ", fluxCalObject
    print "\n---------------------\n"

    obs = ObsFile(obsFileName)
    obs.loadWvlCalFile(wvlCalFileName)
    obs.loadFlatCalFile(flatCalFileName)
    obs.loadFluxCalFile(fluxCalFileName)
    print "loaded data file and calibrations\n---------------------\n"

    HotPixFile = getTimeMaskFileName(obsFileName)
    if not os.path.exists(HotPixFile):
        hp.findHotPixels(obsFileName,HotPixFile)
        print "Flux file pixel mask saved to %s"%(HotPixFile)
    obs.loadHotPixCalFile(HotPixFile)
    print "Hot pixel mask loaded %s"%(HotPixFile)

    nRow = obs.nRow
    nCol = obs.nCol
    obsTime = obs.getFromHeader("exptime")
    wvlBinEdges,obsSpectra,obsIntTime = loadSpectra(obs,nCol,nRow)
    nWvlBins=len(wvlBinEdges)-1
    #divide spectrum by bin widths
    binWidths = np.empty(nWvlBins)
    for i in xrange(nWvlBins):
        binWidths[i] = wvlBinEdges[i+1]-wvlBinEdges[i]
    obsSpectra = obsSpectra/binWidths
    #divide by effective integration time
    for iRow in xrange(nRow):
        for iCol in xrange(nCol):
            obsSpectra[iRow][iCol] = obsSpectra[iRow][iCol] / obsIntTime[iRow][iCol]

    #print np.shape(obsSpectra)
    #print nRow
    #print nCol
    #print nWvlBins

    medianObsSpectrum = calculateMedian(obsSpectra,nCol,nRow,nWvlBins)
    print "target spectrum loaded\n---------------------\n"

    if skyfile != "None":
        sky = ObsFile(skyFileName)
        sky.loadWvlCalFile(wvlCalFileName)
        sky.loadFlatCalFile(flatCalFileName)
        sky.loadFluxCalFile(fluxCalFileName)
        skyTime = sky.getFromHeader("exptime")

        skyHotPixFile = getTimeMaskFileName(skyFileName)
        if not os.path.exists(skyHotPixFile):
            hp.findHotPixels(skyFileName,skyHotPixFile)
            print "Flux file pixel mask saved to %s"%(skyHotPixFile)
        obs.loadHotPixCalFile(HotPixFile)
        print "Hot pixel mask loaded %s"%(skyHotPixFile)

        wvlBinEdges,skySpectra,skyIntTime = loadSpectra(sky,nCol,nRow)

        skySpectra = skySpectra/binWidths
        #divide by effective integration time
        for iRow in xrange(nRow):
            for iCol in xrange(nCol):
                skySpectra[iRow][iCol] = skySpectra[iRow][iCol] / skyIntTime[iRow][iCol]

        skySpectrum = calculateMedian(skySpectra, nCol, nRow, nWvlBins)
        #skySpectrum = skySpectrum*float(obsTime)/float(skyTime) #scale sky spectrum to target observation time
        print "sky spectrum loaded\n---------------------\n"
    else:
        #if no sky file given, estimate sky spectrum as median spectrum of obs file, assuming object is tiny
        skySpectrum = medianObsSpectrum
        print "sky spectrum estimated as median of target file spectrum\n---------------------\n"

    #subtract sky spectrum from every pixel
    allSkySpectrum = obsSpectra-skySpectrum
    #set any negative values to 0 after sky subtraction
    allSkySpectrum[allSkySpectrum<0]=0

    #take median of remaining sky subtracted spectra to get median object spectrum
    finalSpectrum = calculateMedian(allSkySpectrum,nCol,nRow,nWvlBins)
    
    #load std spectrum for comparison
    realSpectra = loadStd(objectName,wvlBinEdges)
    print "real std spectrum loaded for reference\n---------------------\n"

    #create plots
    plotDir = "/home/srmeeker/ARCONS-pipeline/fluxcal/test/plots"
    plotFileName = "%s_from_%s%s.pdf"%(objectName,fluxCalObject,filenum)
    fullFluxPlotFileName = os.path.join(plotDir,plotFileName)
    pp = PdfPages(fullFluxPlotFileName)
    matplotlib.rcParams['font.size']=6

    #calculate midpoints of wvl bins for plotting
    wvls = np.empty((nWvlBins),dtype=float)
    for n in xrange(nWvlBins):
        binsize=wvlBinEdges[n+1]-wvlBinEdges[n]
        wvls[n] = (wvlBinEdges[n]+(binsize/2.0))

    plt.figure()

    ax1 = plt.subplot(231)
    ax1.set_title('ARCONS median flat/flux cal\'d obs in counts')
    plt.plot(wvls,medianObsSpectrum)
    #plt.show()
    ax2 = plt.subplot(232)
    ax2.set_title('ARCONS median flat/flux cal\'d sky in counts')
    plt.plot(wvls,skySpectrum)
    #plt.show()
    ax5 = plt.subplot(233)
    ax5.set_title('Sensitivity Spectrum')
    plt.plot(wvls,obs.fluxWeights)
    ax3 = plt.subplot(234)
    ax3.set_title('MKID data minus sky in counts')
    plt.plot(wvls,finalSpectrum/max(finalSpectrum))
    ax4 = plt.subplot(235)
    ax4.set_title('Rebinned Std Spectrum of %s'%(objectName))
    plt.plot(wvls,realSpectra)

    #ax7 = plt.subplot(337)
    #ax7.set_title('Flux Cal\'d ARCONS Spectrum of Std')
    #plt.plot(wvls,fluxFactors*subtractedSpectra)

    pp.savefig()
    pp.close()

    #del obs
    #del sky

    print "output plots to %s\n---------------------\n"%(fullFluxPlotFileName)

    txtDir = "/home/srmeeker/ARCONS-pipeline/fluxcal/test/txt"
    txtFileName = "%s_from_%s%s.txt"%(objectName,fluxCalObject,filenum)
    fullFluxTxtFileName = os.path.join(txtDir,txtFileName)

    outarr = np.empty((len(medianObsSpectrum),2),dtype=float)
    outarr[:,0]=wvls
    outarr[:,1]=medianObsSpectrum
    #save sensitivity spectrum to file
    np.savetxt(fullFluxTxtFileName, outarr)
    
    print "output txt file to %s\n---------------------\n"%(fullFluxPlotFileName)


if __name__ == '__main__':
    main()














