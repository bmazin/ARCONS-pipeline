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

    #wvldir = "/Scratch/waveCalSolnFiles/oldbox_numbers/20121205"
    #objectName = "crabNight1"

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
    

    #GET RAW PIXEL COUNT IMAGE TO CALCULATE CORRECTION FACTORS
    print "Making raw cube to get dead time correction"
    cubeDict = obs.getSpectralCube(firstSec=startTime, integrationTime=intTime, weighted=False)

    cube= np.array(cubeDict['cube'], dtype=np.double)
    wvlBinEdges= cubeDict['wvlBinEdges']
    effIntTime= cubeDict['effIntTime']
    print "median effective integration time = ", np.median(effIntTime)

    nWvlBins=len(wvlBinEdges)-1
    print "cube shape ", np.shape(cube)
    print "effIntTime shape ", np.shape(effIntTime)

    #add third dimension to effIntTime for  broadcasting
    effIntTime = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
    #put cube into counts/s in each pixel
    cube /= effIntTime

    #CALCULATE DEADTIME CORRECTION
    #NEED TOTAL COUNTS PER SECOND FOR EACH PIXEL TO DO PROPERLY
    #ASSUMES SAME CORRECTION FACTOR APPLIED FOR EACH WAVELENGTH, MEANING NO WL DEPENDANCE ON DEAD TIME EFFECT
    DTCorr = np.zeros((np.shape(cube)[0],np.shape(cube)[1]),dtype=float)
    for f in range(0,np.shape(cube)[2]):
        print cube[:,:,f]
        print '-----------------------'
        DTCorr += cube[:,:,f]
        print DTCorr
        print '\n=====================\n'
    #Correct for 100 us dead time
    DTCorrNew=DTCorr/(1-DTCorr*100e-6)
    CorrFactors = DTCorrNew/DTCorr #This is what the frames need to be multiplied by to get their true values

    print "Dead time correction factors = "
    print CorrFactors

    print "Making Weighted cube"
    #REMAKE CUBE WITH FLAT WEIGHTS AND APPLY DEAD TIME CORRECTION AS WELL

    #load/generate hot pixel mask file
    HotPixFile = getTimeMaskFileName(obsFileName)
    if not os.path.exists(HotPixFile): #check if hot pix file already exists
        hp.findHotPixels(obsFileName,HotPixFile)
        print "Flux file pixel mask saved to %s"%(HotPixFile)
    obs.loadHotPixCalFile(HotPixFile)
    print "Hot pixel mask loaded %s"%(HotPixFile)

    cubeDict = obs.getSpectralCube(firstSec=startTime, integrationTime=intTime, weighted=False)

    cube= np.array(cubeDict['cube'], dtype=np.double)
    wvlBinEdges= cubeDict['wvlBinEdges']
    effIntTime= cubeDict['effIntTime']
    print "median effective integration time = ", np.median(effIntTime)

    nWvlBins=len(wvlBinEdges)-1
    print "cube shape ", np.shape(cube)
    print "effIntTime shape ", np.shape(effIntTime)

    #add third dimension to effIntTime for broadcasting
    effIntTime = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
    #put cube into counts/s in each pixel
    cube /= effIntTime

    #add third dimension to CorrFactors for broadcasting
    CorrFactors = np.reshape(CorrFactors,np.shape(CorrFactors)+(1,))
    #apply dead time correction factors
    cube*=CorrFactors    

    
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
    utils.makeMovie(movieCube,frameTitles=wvls,cbar=True,outName='%s_raw_%s.gif'%(objectName,fileNum), normMin=0, normMax=50)

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














