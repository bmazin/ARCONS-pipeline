from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util.ObsFile import ObsFile
from util.popup import *
from util import MKIDStd
import matplotlib

'''
Demonstration of loading and applying synthetic filters to ObsFile objects
Provides example code for successfully applying filter, and also possible pitfalls

usage:
>>python syntheticFilterTest.py [filter]

options:
[filter] = an optional filter can be provided in command line if you want to play with all
the available options. Default is Johnson V band filter.
'''

#  If you want to make the test plots with a specific filter, provide it in the command line
if len(sys.argv)>1:
    filt = sys.argv[1]
else:
    filt='V'


#  All examples shown here will use data from the W Uma binary, 1SWASP J0002, found here:
fileName = '/ScienceData/PAL2014/20140923/obs_20140924-080031.h5'


#  We'll select a pixel in the middle of the target for our examples
col = 18
row = 16


#  first, load the data into an ObsFile object
obs = ObsFile(fileName)
print "Loaded obsFile", fileName


#  if we tried to simply load a filter now, it would fail as we haven't defined our spectrum's binning.
#  Uncomment the following to see the error you would get.
'''
obs.loadFilter()
'''


#  Wavelength bins can be provided during the filter loading, but the easiest way is to
#  just load the flat calibration first, as that sets everything to the binning used for
#  flat cal. This file already has all its calibrations associated
#  with it in /Scratch/calLookup/lookup.h5. If a file doesn't have it's data written into 
#  this table the following will not work and calibrations must be loaded manually.
obs.loadAllCals()
print "Finished loading calibrations"


#  Now when we load a filter it works fine. A call to loadFilter() with no arguments applies
#  a Johnson V filter by default. To see a list of available filters try calling the function
#  with an unsupported filter. 
#  Uncomment to following to see the error you would get.
'''
obs.loadFilter('this is not a supported filter')
'''


#  You must be careful when using flat/flux calibrations that you don't try to change the
#  spectral binning later. Flat and flux calibrations have a pre-defined binning that the
#  filter must match. That's why it's easiest to just load the calibrations and let the 
#  flat cal binning take over. If you try to override the flat cal binning by loading a 
#  filter with new binning ObsFile will complain when you try to make a spectrum.
#  Uncomment the following to see the error you would get.
'''
myOwnWvlBins = ObsFile.makeWvlBins(energyBinWidth=0.4, wvlStart=4000, wvlStop=11000)
obs.loadFilter(wvlBinEdges = myOwnWvlBins)
obs.getPixelSpectrum(row, col, firstSec=0, integrationTime= -1, weighted=True, fluxWeighted=True)
'''


#  Now that we've covered all the ways the code can fail, let's show examples of it working.
#  We'll load a filter with all of the keywords set to their default. The following command
#  is equivalent to obs.loadFilter().
obs.loadFilter(filterName = 'V', wvlBinEdges = None,switchOnFilter = True)


#  If a different filter is asked for in the command line we'll load it instead
if filt!='V':
    obs.loadFilter(filterName = filt)

#  Once the filter is loaded we can access the transmission curve directly if we want.
#  The filters are loaded from standard curves in MKIDStd, then normalized to 1 and
#  rebinned to whatever wvlBinEdges we provided.
#  Uncomment the following to see a plot of the filter curve and its rebinned shape.
'''
plt.plot(obs.rawFilterWvls, obs.rawFilterTrans, marker='o', markersize = 2, linestyle='-', color='black', label = "Original filter from MKIDStd")
plt.step(obs.filterWvlBinEdges[:-1], obs.filterTrans, where='post', color='blue', label = "Rebinned to wvlBinEdges")
plt.legend()
plt.xlabel('Angstroms')
plt.ylabel('Transmission')
plt.show()
'''


#  By default loadFilter() will also turn on the application of the filter. You can set
#  switchOnFilter to False when calling loadFilter to not activate the filter immediately,
#  and also toggle the filters on and off with switchOffFilter() and switchOnFilter().
obs.switchOffFilter()


#  With the filter off let's get the full spectrum of the selected pixel.
specDict = obs.getPixelSpectrum(row, col, firstSec=0, integrationTime= -1, weighted=True, fluxWeighted=True)
specWvlBinEdges = specDict['wvlBinEdges']
fullSpectrum = specDict['spectrum']


#  Now let's switch the filter back on and overplot the filtered spectrum with the full
#  spectrum.
obs.switchOnFilter()
filteredSpecDict = obs.getPixelSpectrum(row, col, firstSec=0, integrationTime= -1, weighted=True, fluxWeighted=True)
filteredSpecWvlBinEdges = filteredSpecDict['wvlBinEdges']
filteredSpectrum = filteredSpecDict['spectrum']
plt.step(specWvlBinEdges[:-1], fullSpectrum, where = 'post', label = 'Full ARCONS spectrum of single pixel')
plt.step(filteredSpecWvlBinEdges[:-1], filteredSpectrum, where = 'post', label = 'Filtered spectrum')
plt.xlabel('Angstroms')
plt.ylabel('Total counts (not counts/s/Angstrom)')
plt.ylim(0,1200000)
plt.legend()
plt.show()


#  Synthetic filters also apply seamlessly to images. Let's create an image cube that's
#  had a synthetic filter applied and step through the bins to show the filters are 
#  cutting the cube to appropriate levels.
print "Loading full spectral cube for file..."
cubeDict = obs.getSpectralCube(firstSec=0,integrationTime=5,weighted=True,fluxWeighted=True)
cube = np.array(cubeDict['cube'], dtype=np.double)
effIntTime= cubeDict['effIntTime']
cBE = cubeDict['wvlBinEdges']
#add third dimension to effIntTime for broadcasting
effIntTime = np.reshape(effIntTime,np.shape(effIntTime)+(1,))
#put cube into counts/s in each pixel
cube /= effIntTime
print "Scaled spectral cube by effective integration time"
nbins = len(cBE)-1
binCenters = np.zeros((nbins),dtype = int)
for i in xrange(nbins):
    binCenters[i] = cBE[i]+(cBE[i+1]-cBE[i])/2.0
for j in xrange(len(cube[0,0,:])):
    plt.clf()
    plt.step(filteredSpecWvlBinEdges[:-1], filteredSpectrum, where = 'post', label = 'Filtered spectrum')
    plt.axvspan(cBE[j],cBE[j+1], alpha = 0.3, color='green')
    plt.xlabel('Angstroms')
    plt.ylabel('Total counts (not counts/s/Angstrom)')
    plt.ylim(0,1.1*max(filteredSpectrum))
    plt.xlim(3000,13000)
    plt.legend()
    plt.show(block=False)
    form = PopUp(showMe=False,title='%s'%binCenters[j])
    validCube = cube.flatten()[np.isnan(cube.flatten())==False]
    nonZeroCube = validCube[validCube > 0]
    vmax = np.mean(nonZeroCube)+3*np.std(nonZeroCube)
    form.plotArray(cube[:,:,j],vmax = vmax)
    form.show()
    del form

