import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from util.rebin import rebin
import matplotlib
from scipy import interpolate
from scipy.optimize.minpack import curve_fit
from math import exp

"""
Created 2/15/2013 by Seth Meeker
"""
unitdict = {
        'cgs':{ 'h':6.626068e-27,
            'k':1.3806503e-16,
            'c':2.99792458e10},
        'mks':{ 'h':6.626068e-34,
            'k':1.3806503e-23,
            'c':2.99792458e8}
        }
h,k,c = unitdict['cgs']['h'],unitdict['cgs']['k'],unitdict['cgs']['c']

objectName = "HR9087"
#import the known spectrum of the calibrator and rebin to the histogram parameters given
#must be imported into array with dtype float so division later does not have error
std = MKIDStd.MKIDStd()
a = std.load(objectName)
a = std.normalizeFlux(a)
x = a[:,0]
y = a[:,1]

features = [3890,3970,4099,4340,4860,6564,6883,7619]
widths = [50,50,50,50,50,50,50,50]

plt.figure()
ax1 = plt.subplot(1,2,1)
ax1.set_xlim(3000,13000)
plt.plot(x,y)
for feat in features:
   plt.vlines(feat,0,1.5,linestyle='dashed')

#cut out absorption features
#for i in xrange(len(features)):
#    ind = np.where((x<(features[i]+15)) & (x>(features[i]-15)))[0]
#    if len(ind)!=0:
#        ind = ind[len(ind)/2]
#    if y[ind]<y[ind+1] and y[ind]<y[ind-1]:
#        inds = np.where((x >= features[i]+widths[i]) | (x <= features[i]-widths[i]))
#        x = x[inds]
#        y = y[inds]
#plt.plot(x,y,color = 'red')

#set parameters for high wavelength tail interpolation
smooth = 5
fraction = 2.0/3.0
newx = np.arange(int(x[fraction*len(x)]),13000)

#func = interpolate.splrep(x[fraction*len(x):],y[fraction*len(x):],s=smooth)
#func = interpolate.interp1d(x[fraction*len(x):],y[fraction*len(x):],'linear')
#newy = interpolate.splev(newx,func)
#newy = func(newx)

slopeguess = (np.log(y[-1])-np.log(y[fraction*len(x)]))/(x[-1]-x[fraction*len(x)])
print "Guess at exponential slope is %f"%(slopeguess)
guess_a, guess_b, guess_c = float(y[fraction*len(x)]), x[fraction*len(x)], slopeguess
guess = [guess_a, guess_b, guess_c]

fitx = x[fraction*len(x):]
fity = y[fraction*len(x):]

exp_decay = lambda fx, A, x0, t: A * np.exp((fx-x0) * t)
#bb = lambda fx,A,T: 10**11* A * 2*h*c**2 / fx**5 * 1/(np.exp(h*c/(k*T*fx)) - 1)
params, cov = curve_fit(exp_decay, fitx, fity, p0=guess, maxfev=2000)
A, x0, t= params
#A,T = params
print "A = %s\nx0 = %s\nt = %s\n"%(A, x0, t)
best_fit = lambda fx: A * np.exp((fx-x0)*t)
#best_fit = lambda fx:10**11 * A * 2*h*c**2 / fx**5 * 1/(np.exp(h*c/(k*T*fx)) - 1)

calcx = np.array(newx,dtype=float)
newy = best_fit(calcx)

plt.plot(newx,newy,color='green')

finalx = np.concatenate((x,newx[newx>max(x)]))
finaly = np.concatenate((y,newy[newx>max(x)]))

ax2 = plt.subplot(1,2,2)
ax2.set_xlim(3000,13000)

plt.plot(finalx,finaly)
ax1.set_xlim(3000,13000)

for feat in features:
   plt.vlines(feat,0,1.5,linestyle='dashed')

plt.show()



