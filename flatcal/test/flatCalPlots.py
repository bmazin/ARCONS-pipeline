from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
from util.FileName import FileName
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
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
    #obsFileName = FileName(run=run,date='20121210',tstamp='20121211-135628').flat()
    flatCalFileName = FileName(run=run,date='20121210',tstamp='').flatSoln()
    wvlCalFileName = FileName(run=run,date='20121210',tstamp='20121211-133056').calSoln()
    flatCalPath = os.path.dirname(flatCalFileName)
    pp = PdfPages(os.path.join(flatCalPath,'flat20121210.pdf'))
    nPlotsPerPage = 9
    iPlot = 0

    #ob = ObsFile(obsFileName)#('obs_20120919-131142.h5')
    #ob.loadWvlCalFile(wvlCalFileName)#('calsol_20120917-072537.h5')
    #ob.loadFlatCalFile(flatCalFileName)#('flatsol_20120919-131142.h5')
    matplotlib.rcParams['font.size'] = 4

    for testPixelRow in xrange(46):
        for testPixelCol in xrange(44):
    #plot some uncalibrated and calibrated spectra for one pixel
            if sum(cal.spectra[testPixelRow,testPixelCol]) != 0:
                if iPlot % nPlotsPerPage == 0:
                    fig = plt.figure(figsize=(10,10),dpi=100)

                ax = fig.add_subplot(3,3,iPlot%nPlotsPerPage+1)
                ax2 = ax.twinx()

                #flatSpectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,weighted=False)
                #spectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,wvlStart=cal.wvlStart,wvlStop=cal.wvlStop,wvlBinWidth=cal.wvlBinWidth,weighted=False,firstSec=0,integrationTime=-1)

                #weightedSpectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,weighted=True)
                #flatSpectrum,wvlBinEdges = cal.flatFile.getPixelSpectrum(testPixelRow,testPixelCol,wvlStart=cal.wvlStart,wvlStop=cal.wvlStop,wvlBinWidth=cal.wvlBinWidth,weighted=False,firstSec=0,integrationTime=-1)
                flatSpectrum = cal.spectra[testPixelRow,testPixelCol]
                x = cal.wvlBinEdges[0:-1]
                ax.plot(x,cal.wvlMedians,label='median',alpha=.5)
                ax2.plot(x,cal.flatFactors[testPixelRow,testPixelCol,:],'r',label='weights',alpha=.5)
                ax2.set_title('p %d,%d'%(testPixelRow,testPixelCol))
                #ax.plot(x,spectrum+20,label='unweighted spectrum for pixel %d,%d'%(testPixelRow,testPixelCol),alpha=.5)
                #ax.plot(x,weightedSpectrum+10,label='weighted %d,%d'%(testPixelRow,testPixelCol),alpha=.5)
                ax.plot(x,flatSpectrum,label='pixel',alpha=.5)

                ax.legend(loc='lower left')
                ax2.legend(loc='lower right')
                if iPlot%nPlotsPerPage == nPlotsPerPage-1 or (testPixelRow == 45 and testPixelCol == 43):
                    pp.savefig(fig)
                    #plt.show()
                iPlot += 1
    pp.close()

    #display a time-flattened image of the twilight flat as it is and after using itself as it's flat cal
    #cal.flatFile.loadFlatCalFile(flatCalFileName)#('flatsol_20120919-131142.h5')
    #cal.flatFile.displaySec(weighted=True,integrationTime=-1)
    #ob.displaySec(integrationTime=-1)
    #ob.displaySec(weighted=True,integrationTime=-1)


#    for idx in range(0,100,20):
#        factors10 = cal.flatFactors[:,:,idx]
#        plt.matshow(factors10,vmax=np.mean(factors10)+1.5*np.std(factors10))
#        plt.title('Flat weights at %d'%cal.wvlBinEdges[idx])
#        plt.colorbar()
#        plt.savefig('plots/factors%d.png'%idx)
#        plt.show()
if __name__ == '__main__':
    main()
