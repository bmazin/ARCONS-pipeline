import numpy as np
from util.popup import plotArray,PopUp
import os
from photometry.PSFphotometry import PSFphotometry
from util.utils import confirm
from scipy.ndimage.filters import convolve1d


savePath = '/Scratch/dataProcessing/J0337/'
nPhaseBins = 150
wvlStart = 3000 #angstrom
wvlEnd = 6000 #angstrom
#dataPath = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_{:.1f}arcsecAperture.npz'.format(nPhaseBins,wvlStart,wvlEnd,apertureRadius)
dataPath = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_optimalAperture.npz'.format(nPhaseBins,wvlStart,wvlEnd)
print os.path.basename(dataPath)

data = np.load(dataPath)

pixelImages = data['pixelImages']
labels = data['obsTimestamps']

apertures = []
psfDicts = []
for iObs,(img,label) in enumerate(zip(pixelImages,labels)):
    centroidGuess = np.unravel_index(np.argmax(img),np.shape(img))
    centroidGuess = centroidGuess[::-1] #switch from (row,col) to (x,y)
    psfPhot = PSFphotometry(img,centroid=[centroidGuess],verbose=False,showPlot=False)
    psfDict = psfPhot.PSFfit(aper_radius=5)
    params = psfDict['parameters']
    sigma = psfDict['parameters'][4]
    flux = psfDict['flux']
    optimalApertureRadius = 1.6*sigma
    print int(psfDict['flag']),label,'sigma',sigma,'flux',flux,'aperture',optimalApertureRadius
    psfDict2 = data['psfFits'][iObs]
    print int(psfDict2['flag']),label,'sigma',psfDict2['parameters'][4],'flux',psfDict2['flux'],'aperture',1.6*psfDict2['parameters'][4]
    apertures.append(optimalApertureRadius)
    psfDicts.append(psfDict)

    imgPath = os.path.join(os.path.join(savePath,'pics/'),'apert'+label+'.png')
    pop = PopUp(showMe=False)
    pop.plotArray(img,title=label+' \sigma={},f={}'.format(sigma,flux),vmax=np.max(img))
    pop.fig.savefig(imgPath)

apertures = np.array(apertures)


