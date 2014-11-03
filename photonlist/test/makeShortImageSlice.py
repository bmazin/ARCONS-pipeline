import matplotlib.pyplot as mpl
import numpy as np
from util import FileName
from util import utils
import photonlist.photlist as pl
import photonlist.RADecImage as rdi


def makeShortImageSlice(vImage=None,
                        fileName1 = FileName.FileName(run='PAL2012',date='20121211', tstamp='20121212-045902').photonList(),
                        fileName2 = FileName.FileName(run='PAL2012',date='20121211', tstamp='20121212-060934').photonList(),
                        firstSec1=0, integrationTime1=30,
                        firstSec2=0, integrationTime2=30,
                        doWeighted=True, pcliplo=1.0, pcliphi = 98.,
                        wvlMin=3500, wvlMax=12000):
    '''
    Creates a figure used in the ARCONS pipeline paper. Make two versions of a
    coadded pair of short (30s) image slices of the Crab
    pulsar in RA/dec coordinates one weighted by exposure time,
    and one not.
    
    INPUTS:
        vImage - provide return value from a previous run of this function,
                so that it won't have to recalculate everything (note that
                this is *only* for re-displaying. E.g. won't do doWeighted
                differently regardless of how you set it.
        fileName1, fileName2 - two input photon list file names
        firstSec1, firstSec2 - start time (sec) in each of the two photon list files
        integrationTime1, integrationTime2 - integration time to include for each of the two files.
        doWeighted - if True, weight by flat and flux calibrations
        pcliplo, pcliphi - percentile clips to set the color scale stretch
                for the two images, based on the exposure time weighted one.
        wvlMin, wvlMax - minimum and maximum wavelengths to include.
        
    OUTPUTS:
        Plots to screen, and saves the image as 'expTimeWeighting.eps'
    '''
    
    if vImage is None:
        phList1 = pl.PhotList(fileName1)
        phList2 = pl.PhotList(fileName2)
        virtualImage = rdi.RADecImage(phList1,nPixRA=400,nPixDec=400,cenRA=1.4596725441339724,
                                      cenDec=0.38422539085925933,firstSec=firstSec1,
                                      integrationTime=integrationTime1,doWeighted=doWeighted, 
                                      wvlMin=wvlMin, wvlMax=wvlMax)
        virtualImage.loadImage(phList2,firstSec=firstSec2,integrationTime=integrationTime2,
                          doWeighted=doWeighted,wvlMin=wvlMin,wvlMax=wvlMax,
                          doStack=True)
    else:
        virtualImage = vImage
    
    expWtIm = virtualImage.image*virtualImage.expTimeWeights
    vmin = np.percentile(expWtIm[np.isfinite(expWtIm)], pcliplo)
    vmax = np.percentile(expWtIm[np.isfinite(expWtIm)], pcliphi)
    
    
    fig = mpl.figure(figsize=(11,7))
    mpl.subplot(121)
    virtualImage.display(expWeight=False, normMin=vmin, normMax=vmax, colormap=mpl.cm.gray,
                         noAxis=True)
    ax1 = mpl.gca()
    
    mpl.subplot(122)
    virtualImage.display(expWeight=True, normMin=vmin, normMax=vmax, colormap=mpl.cm.gray,
                         noAxis=True)
    ax2 = mpl.gca()

    #cax = fig.add_axes([0.9, 0.1, 0.03, 0.8])
    
    mpl.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.91, wspace=0.05, hspace=0.05)
    cb = mpl.colorbar(ax=[ax1,ax2], orientation='horizontal', pad=0.04, shrink=0.8, 
                      fraction=0.1, aspect=40)
    cb.set_label('Counts')
    
    mpl.savefig('expTimeWeighting.eps', bbox_inches='tight', pad_inches=0)
    
    return virtualImage
