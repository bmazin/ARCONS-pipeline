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
from testImageStack import makeImageStack
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib

"""

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

    objectName = "hz21"
    fileNum=0
    energyBinWidth = 0.1
    #make bins for 3000 to 13000
    wvlStart = 3000
    wvlStop = 13000
    wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth,wvlStart,wvlStop)
    nWvlBins = len(wvlBinEdges)-1
    binWidths = np.empty(nWvlBins)
    print "Showing bin widths for %i bins"%(nWvlBins)
    for i in xrange(nWvlBins):
        binWidths[i] = wvlBinEdges[i+1]-wvlBinEdges[i]
    print binWidths

    nVirtPixX=250
    nVirtPixY=250
    cube = np.zeros((nVirtPixX,nVirtPixY,nWvlBins),dtype=float)

    for n in xrange(nWvlBins):
        print "Making image for wvls %i to %i"%(wvlBinEdges[n], wvlBinEdges[n+1])
        virtualImage, imageStack, medImage = makeImageStack(fileNames='photons_*.h5', dir=os.getenv('INTERM_DIR', default="/Scratch")+'/photonLists/20131209',
                   detImage=False, saveFileName='stackedImage.pkl', wvlMin=wvlBinEdges[n],
                   wvlMax=wvlBinEdges[n+1], doWeighted=True, medCombine=False, vPlateScale=0.2,
                   nPixRA=nVirtPixX,nPixDec=nVirtPixY)
        print virtualImage
        print virtualImage.image
        print np.shape(virtualImage.image)
        cube[:,:,n] = virtualImage.image
    
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

    np.savez('%s_raw_%s.npz'%(objectName,fileNum),stack=movieCube,wvls=wvls)
    utils.makeMovie(movieCube,frameTitles=wvls,cbar=True,outName='%s_pl_raw_%s.gif'%(objectName,fileNum), normMin=0, normMax=1000)

    '''
    #load std spectrum for comparison
    try:
        realSpectra = loadStd(objectName,wvlBinEdges)
        print "real std spectrum loaded for reference\n---------------------\n"
        stdTitle = "Rebinned Std Spectrum of %s"%(objectName)
    except KeyError:
        print "Key Error loading MKIDStd"
        realSpectra = np.ones(nWvlBins)
        stdTitle = "No MKIDStd spectrum available for %s"%(objectName)


    #create plots
    plotDir = "/home/srmeeker/ARCONS-pipeline/fluxcal/test/plots"
    plotFileName = "%s_from_%s%s.pdf"%(objectName,fluxCalObject,filenum)
    fullFluxPlotFileName = os.path.join(plotDir,plotFileName)
    pp = PdfPages(fullFluxPlotFileName)
    matplotlib.rcParams['font.size']=6

    plt.figure()

    ax1 = plt.subplot(221)
    ax1.set_title('ARCONS median flat/flux cal\'d obs in counts')
    ax1.set_xlim((4000,11000))
    ax1.set_ylim((min(medianObsSpectrum[(wvls>4000) & (wvls<11000)]),max(medianObsSpectrum[(wvls>4000) & (wvls<8000)])))
    plt.plot(wvls,medianObsSpectrum)
    #plt.show()
    #ax2 = plt.subplot(232)
    #ax2.set_title('ARCONS median flat/flux cal\'d sky in counts')
    #plt.plot(wvls,skySpectrum)
    #plt.show()
    ax5 = plt.subplot(223)
    ax5.set_title('Sensitivity Spectrum')
    ax5.set_xlim((3000,13000))
    ax5.set_ylim((0,5))
    plt.plot(wvls,obs.fluxWeights)
    #ax3 = plt.subplot(234)
    #ax3.set_title('MKID data minus sky in counts')
    #plt.plot(wvls,finalSpectrum/max(finalSpectrum))
    ax4 = plt.subplot(222)
    ax4.set_title(stdTitle)
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
    '''


if __name__ == '__main__':
    main()














