from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from util.utils import rebin
from matplotlib import rcParams
import matplotlib
from scipy import interpolate
from scipy import integrate
from scipy.optimize.minpack import curve_fit
from scipy.constants import h,c,pi
from math import exp

"""
Written by Seth Meeker, 9/9/2013
Test the rebin utility for flux conservation. Find for Vega, rebinning at our PAL2012 array's median energy resolution 
overestimates integrated flux by 2%. Rebinning at other bin widths underestimates by various amounts, also negligible.
Simply run in command line for diagnostic output and example plot: -> python rebinTest.py
Change objectName below to do for different target in MKIDStd library.
"""

def loadFilter(fname):
    #load filter transmission
    try:
        fdata = np.loadtxt(fname,dtype=float)
    except ValueError:
        fdata = np.loadtxt(fname,dtype=float,delimiter = ',')
    filtx = np.array(fdata[:,0])
    filty = np.array(fdata[:,1])
    filty[filty<0.0] = 0.0
    if max(filty)>1.0:
        filty/=100.0
    if max(filtx<2000):
        filtx*=10.0
    #print (filty)

    if fname == 'Rfilter.txt':
        filty = filty/max(filty)*0.8218 #normalize to our Johnson R filter max Transmission
    if fname == 'Vfilter.txt':
        filty = filty/max(filty)*0.8623 #normalize to our Johnson V filter max transmission

    if (filtx[-1] < filtx[0]): #if array goes high to low, reverse the order to the integration does not get a negative
        filtx = filtx[::-1]
        filty = filty[::-1]

    filtWidth = filtx[-1]-filtx[0]
    print "filter width = ",filtWidth
    filtCorrection = integrate.simps(filty,x=filtx)/filtWidth #calculate correction factor for filter width
    print "filter correction = ", filtCorrection
    return filtx, filty, filtWidth, filtCorrection


objectName = "vega"
std = MKIDStd.MKIDStd()
a = std.load(objectName)
stdWvls = np.array(a[:,0],dtype=float) # wavelengths in angstroms
stdFlux = np.array(a[:,1],dtype=float) # flux in counts/s/cm^2/Angs
stdFlux *= (h*c*1E10/stdWvls) #flux in J/s/cm^2/Angs or W/cm^2/Angs


#cut spectrum arrays to our wavelength region
stdFlux = stdFlux[(stdWvls>3000) & (stdWvls<13000)]
stdWvls = stdWvls[(stdWvls>3000) & (stdWvls<13000)]


plt.plot(stdWvls[(stdWvls>3158)&(stdWvls<11089)],stdFlux[(stdWvls>3158)&(stdWvls<11089)])
total = integrate.simps(stdFlux[(stdWvls>3158)&(stdWvls<11089)],x=stdWvls[(stdWvls>3158)&(stdWvls<11089)])
print "Original total flux = ", total
#rebin test
dE = 0.3936 #eV
start = 3000 #Angs
stop = 13000 #Angs
for n in xrange(10):
    n=float(n)+1.0
    enBins = ObsFile.makeWvlBins(dE/n,start,stop)
    rebinned = rebin(stdWvls,stdFlux,enBins)
    re_wl = rebinned[:,0]
    re_flux = rebinned[:,1]
    print re_wl[0]
    print re_wl[-1]
    plt.plot(re_wl[(re_wl>3158)&(re_wl<11089)],re_flux[(re_wl>3158)&(re_wl<11089)])
    total = integrate.simps(re_flux[(re_wl>3158)&(re_wl<11089)],x=re_wl[(re_wl>3158)&(re_wl<11089)])
    print "Total flux for dE = ",dE/n, "eV is ", total , "W/cm^2"
plt.show()

