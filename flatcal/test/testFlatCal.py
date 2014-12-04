from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
from util.FileName import FileName
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages

def main():
    np.set_printoptions(threshold=np.nan)
    testPixelRow = 24
    testPixelCol = 17
    #obs_20120919-131142.h5,obs_20120919-131346.h5
    #create a cal file from a twilight flat
    cal = FlatCal('../../params/flatCal.dict')
    #open another twilight flat as an observation and apply a wavelength cal and the new flat cal
#    run='LICK2012'
#    obsFileName = FileName(run=run,date='20120918',tstamp='20120919-131142').flat()
#    flatCalFileName = FileName(run=run,date='20120918',tstamp='20120919-131448').flatSoln()
#    wvlCalFileName = FileName(run=run,date='20120916',tstamp='20120917-072537').calSoln()

    run = 'PAL2012'
    obsFileName = FileName(run=run,date='20121211',tstamp='20121212-140003').obs()
    flatCalFileName = FileName(run=run,date='20121210',tstamp='').flatSoln()
    wvlCalFileName = FileName(run=run,date='20121210',tstamp='20121211-133056').calSoln()
    flatCalPath = os.path.dirname(flatCalFileName)

    ob = ObsFile(obsFileName)#('obs_20120919-131142.h5')
    ob.loadWvlCalFile(wvlCalFileName)#('calsol_20120917-072537.h5')
    ob.loadFlatCalFile(flatCalFileName)#('flatsol_20120919-131142.h5')

    #plot some uncalibrated and calibrated spectra for one pixel
    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    print ob.getPixelCount(testPixelRow,testPixelCol)

    #flatSpectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,weighted=False)
    spectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,wvlStart=cal.wvlStart,wvlStop=cal.wvlStop,wvlBinWidth=cal.wvlBinWidth,weighted=False,firstSec=0,integrationTime=-1)

    weightedSpectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,weighted=True)
    #flatSpectrum,wvlBinEdges = cal.flatFile.getPixelSpectrum(testPixelRow,testPixelCol,wvlStart=cal.wvlStart,wvlStop=cal.wvlStop,wvlBinWidth=cal.wvlBinWidth,weighted=False,firstSec=0,integrationTime=-1)
    flatSpectrum = cal.spectra[testPixelRow,testPixelCol]
    x = wvlBinEdges[0:-1]
    ax.plot(x,cal.wvlMedians,label='median spectrum',alpha=.5)
    ax2.plot(x,cal.flatFactors[testPixelRow,testPixelCol,:],label='pixel weights',alpha=.5)
    ax2.set_title('flat weights for pixel %d,%d'%(testPixelRow,testPixelCol))
    ax.plot(x,spectrum+20,label='unweighted spectrum for pixel %d,%d'%(testPixelRow,testPixelCol),alpha=.5)
    ax.plot(x,weightedSpectrum+10,label='weighted %d,%d'%(testPixelRow,testPixelCol),alpha=.5)
    ax.plot(x,flatSpectrum+30,label='flatFile %d,%d'%(testPixelRow,testPixelCol),alpha=.5)

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3),fancybox=True,ncol=3)
    plt.show()

    #display a time-flattened image of the twilight flat as it is and after using itself as it's flat cal
    #cal.flatFile.loadFlatCalFile(flatCalFileName)#('flatsol_20120919-131142.h5')
    #cal.flatFile.displaySec(weighted=True,integrationTime=-1)
    #ob.displaySec(integrationTime=-1)
    #ob.displaySec(weighted=True,integrationTime=-1)


    for idx in range(0,100,20):
        factors10 = cal.flatFactors[:,:,idx]
        plt.matshow(factors10,vmax=np.mean(factors10)+1.5*np.std(factors10))
        plt.title('Flat weights at %d'%cal.wvlBinEdges[idx])
        plt.colorbar()
        plt.savefig('plots/factors%d.png'%idx)
        plt.show()
if __name__ == '__main__':
    main()
