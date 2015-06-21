import numpy as np
from util.popup import plotArray,PopUp
import os
from photometry.PSFphotometry import PSFphotometry
from util.utils import confirm
from scipy.ndimage.filters import convolve1d


savePath = '/Scratch/dataProcessing/J0337/'
nPhaseBins = 150
wvlStart = 3000 #angstrom
wvlEnd = 8000 #angstrom
apertureRadius=3.#arcsec
dataPath = '/Scratch/dataProcessing/J0337/profiles2014_{}bins_{}-{}angstroms_{:.1f}arcsecAperture.npz'.format(nPhaseBins,wvlStart,wvlEnd,apertureRadius)
print os.path.basename(dataPath)

data = np.load(dataPath)

pixelImages = data['pixelImages']
labels = data['obsTimestamps']

apertures = []
psfDicts = []
for img,label in zip(pixelImages,labels):
    centroidGuess = np.unravel_index(np.argmax(img),np.shape(img))
    centroidGuess = centroidGuess[::-1] #switch from (row,col) to (x,y)
    psfPhot = PSFphotometry(img,centroid=[centroidGuess],verbose=False,showPlot=False)
    psfDict = psfPhot.PSFfit(aper_radius=5)
    params = psfDict['parameters']
    sigma = psfDict['parameters'][4]
    flux = psfDict['flux']
    optimalApertureRadius = 1.6*sigma
    print int(psfDict['flag']),label,'sigma',sigma,'flux',flux,'aperture',optimalApertureRadius
    apertures.append(optimalApertureRadius)
    psfDicts.append(psfDict)

#    imgPath = os.path.join(os.path.join(savePath,'pics/'),'apert'+label+'.png')
#    pop = PopUp(showMe=False)
#    pop.plotArray(img,title=label+' \sigma={},f={}'.format(sigma,flux),vmax=np.max(img))
#    pop.fig.savefig(imgPath)

apertures = np.array(apertures)

def smooth(x, nPts=8.):
    box = 1.*np.ones(nPts)/nPts
    return convolve1d(x,box)

dates = ['20140925','20140926','20141021']
listsByDate = []
smoothedLists = []

for date in dates:
    dateMask = np.array([item.startswith(date) for item in labels])
    apertsOnDate = apertures[dateMask]
    smoothAperts = smooth(apertsOnDate)

    listsByDate.append(apertsOnDate)
    smoothedLists.append(smoothAperts)

smoothedAperts = np.concatenate(smoothedLists)

dataPathApert = '/Scratch/dataProcessing/J0337/apertures2014_{}bins_{}-{}angstroms_{:.1f}arcsecAperture_sigma.npz'.format(nPhaseBins,wvlStart,wvlEnd,apertureRadius)
apertures = np.array(apertures)

np.savez(dataPathApert,optimalApertureRadii=apertures,smoothOptimalApertureRadii=smoothedAperts,psfDicts=psfDicts,**data)

aperturePath = 'apertureRadiusList.txt'
with open(aperturePath,'w') as apertureFile:
    for label,(apert,smoothApert) in zip(labels,zip(apertures,smoothedAperts)):
        apertureFile.write('{}\t{}\t{}\n'.format(label,apert,smoothApert))
            
