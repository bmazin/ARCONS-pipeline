from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
import matplotlib.pyplot as plt
import numpy as np

def main():
    np.set_printoptions(threshold=np.nan)
    #obs_20120919-131142.h5,obs_20120919-131346.h5
    cal = FlatCal('obs_20120919-131346.h5','calsol_20120917-072537.h5','flatsol_20120919-131142.h5')
    ob = ObsFile('obs_20120919-131142.h5')
    ob.loadWvlCalFile('calsol_20120917-072537.h5')
    ob.loadFlatCalFile('flatsol_20120919-131142.h5')

    fig = plt.figure()
    ax = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)
    testPixelRow = 17
    testPixelCol = 12
    print ob.getPixelCount(testPixelRow,testPixelCol)

    #flatSpectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,weighted=False)
    spectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,wvlStart=cal.wvlStart,wvlStop=cal.wvlStop,wvlBinWidth=cal.wvlBinWidth,weighted=False,firstSec=0,integrationTime=-1)

    weightedSpectrum,wvlBinEdges = ob.getPixelSpectrum(testPixelRow,testPixelCol,weighted=True)
    flatSpectrum,wvlBinEdges = cal.flatFile.getPixelSpectrum(testPixelRow,testPixelCol)
    x = wvlBinEdges[0:-1]
    ax.plot(x,cal.wvlMedians,label='median spectrum')
    ax2.plot(x,cal.flatFactors[testPixelRow,testPixelCol,:],label='pixel weights')
    ax.plot(x,spectrum,label='unweighted %d,%d'%(testPixelRow,testPixelCol))
    ax.plot(x,weightedSpectrum,label='weighted %d,%d'%(testPixelRow,testPixelCol))
    ax.plot(x,flatSpectrum,label='flatFile %d,%d'%(testPixelRow,testPixelCol))

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.3),fancybox=True,ncol=3)
    plt.show()
    cal.flatFile.loadFlatCalFile('flatsol_20120919-131142.h5')
    cal.flatFile.displaySec(weighted=True,integrationTime=-1)
    ob.displaySec(integrationTime=-1)
    ob.displaySec(weighted=True,integrationTime=-1)


if __name__ == '__main__':
    main()
