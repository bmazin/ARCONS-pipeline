import numpy as np
import tables
import matplotlib.pyplot as plt
from util.FileName import FileName
from util.popup import plotArray,PopUp
import multiprocessing
from photonlist.photlist import writePhotonList,PhotList,xyPack
from headers.ArconsHeaders import PulsarPhotonList
from timing import photonTiming
import ephem
import functools
from util.readDict import readDict
from util.ObsFile import ObsFile

RADIANS_PER_ARCSEC = np.pi/648000.

def extractAperturePhotons(plPath,centroidRaStr='03:37:43.826',centroidDecStr='14:15:14.828',apertureRadius=1.,nPhaseBins=150,wvlRange=(3000,13000)):
    print 'reading',plPath
    centroidRa = float(ephem.hours(centroidRaStr))
    centroidDec = float(ephem.degrees(centroidDecStr))
    pl = PhotList(plPath)
    apertureRadius *= RADIANS_PER_ARCSEC
    wvlStart,wvlEnd = wvlRange
    print 'reading photons'
    photons = pl.photTable.readWhere('(sqrt((ra-centroidRa)**2. + (dec-centroidDec)**2.) < apertureRadius)')
    del pl
    print len(photons),'read'
    phases = photons['pulsePhase']
    wavelengths = photons['wavelength']
    print 'done reading'
    del photons
    return {'phases':phases,'wavelengths':wavelengths}


if __name__=='__main__':
    savePath = '/Scratch/dataProcessing/J0337/'
    paramFile = 'j0337.dict'
    params = readDict()
    params.read_from_file(paramFile)
    run = params['run']
    sunsetDates = []
    flatDates = []
    obsSequences = []
    parFiles = []

    sunsetDates.append(params['sunsetDate0'])
    flatDates.append(params['flatDate0'])
    parFiles.append(params['parFile0'])
    obsSequences.append(params['obsSequence0'])
    
    sunsetDates.append(params['sunsetDate1'])
    flatDates.append(params['flatDate1'])
    parFiles.append(params['parFile1'])
    obsSequences.append(params['obsSequence1'])

    sunsetDates.append(params['sunsetDate2'])
    flatDates.append(params['flatDate2'])
    parFiles.append(params['parFile2'])
    obsSequences.append(params['obsSequence2'])

    plPaths = []
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        sunsetDate = sunsetDates[iSeq]
        plPaths.append([FileName(run=run,date=sunsetDate,tstamp=ts).photonList() for ts in obsSequence])

    obsTimestamps = np.concatenate(obsSequences)
    plPaths = [path for pathList in plPaths for path in pathList]

    nPhaseBins = 1500 #bin it finely, we can rebin later
    wvlStart = 3000 #angstrom
    wvlEnd = 13000 #angstrom
    wvlRange = (wvlStart,wvlEnd)
    centroidRaStr='03:37:43.826'
    centroidDecStr='14:15:14.828'


#    #open up the first photon list and extract the wavelength bin edges used to make the flatcal and fluxcal
#    pl0 = PhotList(plPaths[0])
#    wvlBinEdges = np.array(pl0.file.root.flatcal.wavelengthBins.read())
#    del pl0


#    apertureData = np.loadtxt('smoothApertList.txt',dtype={'names':('obsTimestamp','oldApert','apertureRadius'),'formats':('S80','f8','f8')},delimiter='\t')
#    apertureRadiusDict = {eachItem['obsTimestamp']:eachItem['apertureRadius'] for eachItem in apertureData}
#    apertureList = [apertureRadiusDict[timestamp] for timestamp in obsTimestamps]
    

    #reduce arguments to just plPath and apertureRadius
    extractAperturePhotonsTwoArgs = functools.partial(extractAperturePhotons,nPhaseBins=nPhaseBins,wvlRange=wvlRange,centroidRaStr=centroidRaStr,centroidDecStr=centroidDecStr)
    #combine these two arguments into one tuple
    def extractAperturePhotonsOneArg(foldArg):
        plPath,apertureRadius = foldArg
        return extractAperturePhotonsTwoArgs(plPath,apertureRadius=apertureRadius)

    dataPathApert = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_{:.1f}arcsecAperture_sigma.npz'.format(150,3000,8000,3.)
    apertureDataDict = np.load(dataPathApert)
    apertureList = apertureDataDict['smoothOptimalApertureRadii']
    tstamps = apertureDataDict['obsTimestamps']
    psfDicts = apertureDataDict['psfDicts']
    tstampMask = np.array([tstamp in obsTimestamps for tstamp in tstamps])
    apertureList = apertureList[tstampMask]
    tstamps = tstamps[tstampMask]
    psfDicts = psfDicts[tstampMask]

    extractFuncArgs = zip(plPaths,apertureList)

    pool = multiprocessing.Pool(processes=1)
    outDicts = pool.map(extractAperturePhotonsOneArg,extractFuncArgs)

    phases = np.concatenate([d['phases'] for d in outDicts])
    wavelengths = np.concatenate([d['wavelengths'] for d in outDicts])
    
    print 'sorting phases'
    sortedIndices = np.argsort(phases)
    phases = phases[sortedIndices]
    wavelengths = wavelengths[sortedIndices]
    nPhotons = len(phases)
    print len(phases),'photons found'

    photonRecArray = np.recarray(dtype=[('phase',np.double),('wavelength',np.double)],shape=(nPhotons,))
    photonRecArray['phase'] = phases
    photonRecArray['wavelength']=wavelengths

    del phases,wavelengths

    fullPhotonListFileName = '/Scratch/dataProcessing/J0337/masterPhotons3.h5'
    plFile = tables.openFile(fullPhotonListFileName, mode='w')
    plGroup = plFile.createGroup('/', 'photons', 'Group containing photon list')
    #make a description for the table we will put in photFile
    desc = {'phase':tables.Float64Col(),'wavelength':tables.Float64Col()}
    zlibFilter = tables.Filters(complevel=2, complib='zlib', fletcher32=False)
    photTable = plFile.createTable(plGroup,'photTable',title='Phase,Wavelength pairs',description=desc,expectedrows=nPhotons,filters=zlibFilter)
    print 'append'
    photTable.append(photonRecArray)
    #photTable.cols.wavelength = wavelengths[:]
    #photTable.cols.phase = phases[:]
    print 'index wavelength'
    photTable.cols.wavelength.createCSIndex()

    plFile.close()
    
