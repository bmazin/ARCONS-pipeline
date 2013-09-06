import numpy as np
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

plt.plot(wvls,curveThru)
plt.xlim(3900,13000)
plt.ylim(0,0.1)

outarr = np.empty((len(wvls),2),dtype=float)
outarr[:,0]=wvls
outarr[:,1]=curveThru
#save throughput curve to file
np.savetxt("throughput.txt", outarr)

peak = max(maxThru[(wvls>3900) & (wvls<6000)])
peakwvl = wvls[maxThru==max(maxThru[(wvls>3900) & (wvls<6000)])]

curveThru /= curveThru[wvls == peakwvl]
#print curveThru

curveThru *= peak
#print curveThru

plt.plot(wvls,curveThru)
plt.xlim(3900,13000)
plt.ylim(0,0.1)
plt.show()

outarr = np.empty((len(wvls),2),dtype=float)
outarr[:,0]=wvls
outarr[:,1]=curveThru
#save throughput curve to file
np.savetxt("throughput_high.txt", outarr)
