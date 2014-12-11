import numpy as np
import tables
import matplotlib.pyplot as plt
from util.FileName import FileName
from util.popup import plotArray,PopUp
import multiprocessing
from photonlist.photlist import writePhotonList,PhotList
from headers.ArconsHeaders import PulsarPhotonList
from timing import photonTiming
import ephem

RADIANS_PER_ARCSEC = np.pi/648000.

def plotPulseProfile(ax,phaseBinEdges,pulseProfile,**kwargs):
    doublePhaseBinEdges = np.append(phaseBinEdges,phaseBinEdges[1:]+1.)
    doublePulseProfile = np.append(pulseProfile,pulseProfile)
    ax.plot(doublePhaseBinEdges,np.append(doublePulseProfile,doublePulseProfile[-1]),drawstyle='steps-post',**kwargs)

def foldPL(plPath,centroidRaStr='03:37:43.826',centroidDecStr='14:15:14.828',apertureRadius=1.,nPhaseBins=160):
    print 'folding',plPath
    centroidRa = float(ephem.hours(centroidRaStr))
    centroidDec = float(ephem.degrees(centroidDecStr))
    pl = PhotList(plPath)
    apertureRadius *= RADIANS_PER_ARCSEC
    print 'reading photons'
    photons = pl.photTable.readWhere('sqrt((ra-centroidRa)**2. + (dec-centroidDec)**2.) < apertureRadius')
    
    print 'binning phases'
    del pl

    phases = photons['pulsePhase']
    print phases[0:10]
    phaseProfile,phaseBinEdges = np.histogram(phases,bins=nPhaseBins,range=(0.,1.))
    print phaseProfile
    profileErrors = np.sqrt(phaseProfile)
    print 'done folding'
    del photons,phases
    return {'phaseProfile':phaseProfile,'phaseBinEdges':phaseBinEdges,'profileErrors':profileErrors}


if __name__=='__main__':
    savePath = '/Scratch/dataProcessing/J0337/'
    run = 'PAL2014'
    sunsetDates = []
    flatDates = []
    obsSequences = []
    parFiles = []
    

    obsSequence0 = [
    '20140925-111532',
    '20140925-112035',
    '20140925-112538',
    '20140925-113041',
    '20140925-113544',
    '20140925-114047',
    '20140925-114550',
    '20140925-115053',
    '20140925-115605',
    '20140925-120112',
    '20140925-120615',
    '20140925-121118',
    '20140925-121621',
    '20140925-122124',
    '20140925-122627',
    '20140925-123130',
    '20140925-123633'
    ]
    sunsetDates.append('20140924')
    flatDates.append('20140924')
    obsSequences.append(obsSequence0)
    parFiles.append('J0337-mjd-56924.par')
    

    obsSequence1 = [
    '20140926-081832',
    '20140926-082336',
    '20140926-082840',
    '20140926-083344',
    '20140926-083848',
    '20140926-084352',
    '20140926-084856',
    '20140926-085400',
    '20140926-085904',
    '20140926-090408',
    '20140926-090912',
    '20140926-091416',
    '20140926-091920',
    '20140926-092424',
    '20140926-092928',
    '20140926-093432',
    '20140926-093936',
    '20140926-094440',
    '20140926-094944',
    '20140926-095448',
    '20140926-095952',
    '20140926-100456',
    '20140926-101000',
    '20140926-101504',
    '20140926-102008',
    '20140926-102512',
    '20140926-103016',
    '20140926-103520',
    '20140926-104024',
    '20140926-104528'
    ]
    sunsetDates.append('20140925')
    flatDates.append('20140924')
    obsSequences.append(obsSequence1)
    parFiles.append('J0337-mjd-56924.par')

    obsSequence2 = [
    '20141021-063246',
    '20141021-063750',
    '20141021-064253',
    '20141021-064758',
    '20141021-065303',
    '20141021-065807',
    '20141021-070311',
    '20141021-071138',
    '20141021-071656',
    '20141021-072200',
    '20141021-072704',
    '20141021-073208',
    '20141021-073713',
    '20141021-074216',
    '20141021-074720',
    '20141021-075225',
    '20141021-075729',
    '20141021-080232',
    '20141021-080737',
    '20141021-081242',
    '20141021-081747',
    '20141021-082252',
    '20141021-082755',
    '20141021-083300',
    '20141021-083805',
    '20141021-084310',
    '20141021-084814',
    '20141021-085318',
    '20141021-085822',
    '20141021-090326',
    '20141021-090832',
    '20141021-091336',
    '20141021-091840',
    '20141021-092345',
    '20141021-092851'
    ]
    sunsetDates.append('20141020')
    flatDates.append('20141020')
    obsSequences.append(obsSequence2)
    parFiles.append('J0337-mjd-56951.par')


    plPaths = []
    for iSeq in range(len(obsSequences)):
        obsSequence = obsSequences[iSeq]
        sunsetDate = sunsetDates[iSeq]
        plPaths.append([FileName(run=run,date=sunsetDate,tstamp=ts).photonList() for ts in obsSequence])

    obsTimestamps = np.concatenate(obsSequences)
    plPaths = [path for pathList in plPaths for path in pathList]

    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()//2)
    outDicts = pool.map(foldPL,plPaths)
    phaseProfiles = np.array([d['phaseProfile'] for d in outDicts])
    profileErrors = np.array([d['profileErrors'] for d in outDicts])
    phaseBinEdges = outDicts[0]['phaseBinEdges']
    savePath = '/Scratch/dataProcessing/J0337/profiles2014.npz'
    np.savez(savePath,phaseProfiles=phaseProfiles,profileErrors=profileErrors,phaseBinEdges=phaseBinEdges,obsTimestamps=obsTimestamps)
    
    fig,ax = plt.subplots()
    plotPulseProfile(ax,phaseBinEdges,phaseProfiles[0])
    print 'done'
    plt.show()
    
