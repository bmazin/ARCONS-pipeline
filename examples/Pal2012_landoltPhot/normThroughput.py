import numpy as np
import sys, os
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from util.rebin import rebin
import matplotlib
from scipy import interpolate
from scipy.optimize.minpack import curve_fit
from math import exp
import matplotlib.pyplot as plt
from matplotlib import cm
from headers import pipelineFlags
import tables

#normalizes a throughput file with the correct curve shape to the correct absolute height as measured by 
#a file that was generated on a better night of data.

maxfile = '/home/srmeeker/scratch/standards/sdss j0926_throughput.npz'

curvefile = '/home/srmeeker/scratch/standards/G158-100_throughput.npz'

maxDict = np.load(maxfile)
curveDict = np.load(curvefile)

maxThru = maxDict['throughput']
curveThru = curveDict['throughput']
wvls = maxDict['wvls']

print wvls
print maxThru
print curveThru

#interpolate over bad data point int throughput
IWvls = np.delete(wvls,-1)
IWvls = np.delete(IWvls,-1)
print IWvls
IcurveThru = np.delete(curveThru,-1)
IcurveThru = np.delete(IcurveThru,-1)
print IcurveThru

'''
exp_decay = lambda fx, A, x0, t: A * np.exp((fx-x0) * t)
slopeguess = (np.log(IcurveThru[-1])-np.log(IcurveThru[-5]))/(IWvls[-1]-IWvls[-5])
guess_a, guess_b, guess_c = float(IcurveThru[-5]), IWvls[-5], slopeguess
guess = [guess_a, guess_b, guess_c]
params, cov = curve_fit(exp_decay, IWvls[-5:], IcurveThru[-5:], p0=guess, maxfev=2000)
A, x0, t= params
print "A = %s\nx0 = %s\nt = %s\n"%(A, x0, t)
best_fit = lambda fx: A * np.exp((fx-x0)*t)
'''

linear = lambda fx, m, x0, b: m * (fx-x0) + b
guess_m, guess_x0, guess_b = (IcurveThru[-2]-IcurveThru[-1])/(IWvls[-2]-IWvls[-1]), 0.0, IcurveThru[-2]
guess = [guess_m, guess_x0, guess_b]
params, cov = curve_fit(linear, IWvls[-5:], IcurveThru[-5:], p0=guess, maxfev=2000)
m, x0, b= params
print "m = %s\nx0 = %s\nb = %s\n"%(m, x0, b)
best_fit = lambda fx: m * (fx-x0) + b


curveThruEnd = best_fit(wvls)[-2:]
print curveThruEnd
curveThru = np.append(IcurveThru, curveThruEnd)

#f = interpolate.interp1d(IWvls,IcurveThru,'linear')
#curveThru = f(wvls)
#curveThru = interpolate.griddata(IWvls,IcurveThru,wvls,'linear')
print "interpolation done"

#plt.plot(wvls,curveThru)
#plt.xlim(3900,13000)
#plt.ylim(0,0.1)

outarr = np.empty((len(wvls),2),dtype=float)
outarr[:,0]=wvls
outarr[:,1]=curveThru
#save throughput curve to file
np.savetxt("throughput.txt", outarr)

#save absolute throughput measurement to fluxcal h5 file
energyBinWidth = 0.1
wvlStart = 3000
wvlStop = 13000
wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth,wvlStart,wvlStop)

fluxFactors = 1.0/curveThru
print "flux Factors = "
print fluxFactors

fluxFlags = np.empty(np.shape(fluxFactors),dtype='int')
fluxFlags.fill(pipelineFlags.fluxCal['good'])   #Initialise flag array filled with 'good' flags. JvE 5/1/2013.
#set factors that will cause trouble to 1
#self.fluxFlags[self.fluxFactors == np.inf] = 1
fluxFlags[fluxFactors == np.inf] = pipelineFlags.fluxCal['infWeight']   #Modified to use flag dictionary - JvE 5/1/2013
fluxFactors[fluxFactors == np.inf]=1.0
fluxFactors[np.isnan(fluxFactors)]=1.0
fluxFlags[np.isnan(fluxFactors)] = pipelineFlags.fluxCal['nanWeight']   #Modified to use flag dictionary - JvE 5/1/2013
fluxFlags[fluxFactors <= 0]=pipelineFlags.fluxCal['LEzeroWeight']   #Modified to use flag dictionary - JvE 5/1/2013
fluxFactors[fluxFactors <= 0]=1.0

fluxCalFileName = "fluxsol_absolute.h5"
if os.path.isabs(fluxCalFileName) == True:
    fullFluxCalFileName = fluxCalFileName
else:
    scratchDir = os.getenv('INTERM_PATH')
    fluxDir = os.path.join(scratchDir,'fluxCalSolnFiles')
    fullFluxCalFileName = os.path.join(fluxDir,fluxCalFileName)

try:
    fluxCalFile = tables.openFile(fullFluxCalFileName,mode='w')
except:
    print 'Error: Couldn\'t create flux cal file, ',fullFluxCalFileName

calgroup = fluxCalFile.createGroup(fluxCalFile.root,'fluxcal','Table of flux calibration weights by wavelength')
caltable = tables.Array(calgroup,'weights',object=fluxFactors,title='Flux calibration Weights indexed by wavelengthBin')
flagtable = tables.Array(calgroup,'flags',object=fluxFlags,title='Flux cal flags indexed by wavelengthBin. 0 is Good')
bintable = tables.Array(calgroup,'wavelengthBins',object=wvlBinEdges,title='Wavelength bin edges corresponding to third dimension of weights array')
fluxCalFile.flush()
fluxCalFile.close()
print "Finished Flux Cal, written to %s"%(fullFluxCalFileName)

plt.plot(wvls,1.0/curveThru)
plt.xlim(3000,13000)
plt.ylim(0,500)
plt.show()

peak = max(maxThru[(wvls>3900) & (wvls<6000)])
peakwvl = wvls[maxThru==max(maxThru[(wvls>3900) & (wvls<6000)])]

curveThru /= curveThru[wvls == peakwvl]
#print curveThru

curveThru *= peak
#print curveThru

outarr = np.empty((len(wvls),2),dtype=float)
outarr[:,0]=wvls
outarr[:,1]=curveThru
#save throughput curve to file
np.savetxt("throughput_high.txt", outarr)
