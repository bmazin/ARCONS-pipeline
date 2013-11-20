from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from util.rebin import rebin
from matplotlib import rcParams
import matplotlib
from scipy import interpolate
from scipy import integrate
from scipy.optimize.minpack import curve_fit
from scipy.constants import h,c,pi
from math import exp

"""
Written by Seth Meeker, 9/9/2013
Calculates telescope throughput at Palomar 200" coude focus
Uses calibrated photodiode (UDT Instruments Silicon Sensor Model 221)
Compare measured power against expected power from standard star (most likely Vega, v=0)
Calculate expected power by multiplying standard spectrum by our V filter transmission,
and by the photodiode absolute response function.
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
    if max(filty)>1.0: #if transmission is in percent, convert to decimal
        filty/=100.0
    if max(filtx<2000): #if wvl is in nm, convert to Angs
        filtx*=10.0
    #print (filty)

    if fname == 'Rfilter.txt':
        filty = filty/max(filty)*0.8218 #normalize to our Johnson R filter max Transmission
    if fname == 'Vfilter.txt':
        filty = filty/max(filty)*0.8623 #normalize to our Johnson V filter max transmission

    if (filtx[-1] < filtx[0]): #if array goes high to low, reverse the order so the integration does not get a negative
        filtx = filtx[::-1]
        filty = filty[::-1]

    filtWidth = filtx[-1]-filtx[0]
    print "filter width = ",filtWidth
    filtCorrection = integrate.simps(filty,x=filtx)/filtWidth #calculate correction factor for filter width
    print "filter correction = ", filtCorrection
    return filtx, filty, filtWidth, filtCorrection


#give measured optometer Amps as argument
if len(sys.argv)<3:
    print "Provide measured Amps from optometer as argument, and desired filter"
    print "Example: >python PalCoudeThroughput.py 1E-4 v"
    print "Available filters: B, V, R
    sys.exit()
AmpsMeas = float(sys.argv[1])
filt = str(sys.argv[2])

#load standard spectrum in counts/s/cm^2/A and convert to watts/cm^2/A
objectName = "vega"
std = MKIDStd.MKIDStd()
a = std.load(objectName)
stdWvls = np.array(a[:,0],dtype=float) # wavelengths in angstroms
stdFlux = np.array(a[:,1],dtype=float) # flux in counts/s/cm^2/Angs
stdFlux *= (h*c*1E10/stdWvls) #flux in J/s/cm^2/Angs or W/cm^2/Angs


#cut spectrum arrays to our wavelength region
stdFlux = stdFlux[(stdWvls>3000) & (stdWvls<13000)]
stdWvls = stdWvls[(stdWvls>3000) & (stdWvls<13000)]


#------test plots------------#
ax1 = plt.subplot(1,1,1)
ax1.set_xlabel("Wavelength in $\AA$")
ax1.set_ylabel("Flux in W/cm$^{2}$/$\AA$")
plt.plot(stdWvls,stdFlux)
plt.show()
#----------------------------#


#multiply by area of Palomar 200" pupil to get watts/Angstrom
# From Sivaramakrishnan et al 2001, D = 5.08m, Ds/D = 0.33
Diam = 508 #cm
Ds = 0.33 * Diam
radp = Diam/2.0
rads = Ds/2.0
area = pi*(radp**2 - rads**2)
print "Primary diameter = ", Diam, " cm"
print "Secondary diameter = ", Ds, " cm"
print "Total collecting area = ", area, " cm^2"
stdFlux*=area #flux in W/Angs

#multiply standard spectrum by V-filter transmission
if filt == ("b" or "B"):
    fname = 'AstrodonB.txt'
if filt == ("v" or "V"):
    fname = 'AstrodonV.txt'
if filt == ("r" or "R"):
    fname = 'AstrodonR.txt'

filtWvls, filtTrans, filtWidth, filtCorrection = loadFilter(fname) #load filter file
filtTransInterp = interpolate.griddata(filtWvls,filtTrans,stdWvls,"linear",fill_value=0)#interpolate filter curve to spectrum wvl spacing
filteredFlux = stdFlux*filtTransInterp

#check total filtered flux matches a v=0 star
checkFlux = filteredFlux/area #W/Angs/cm^2
checkFlux *= 1.0E7 #ergs/s/Angs/cm^2
totalCheckFlux = integrate.simps(checkFlux,x=stdWvls)/filtWidth #ergs/s/cm^2/Angs Vband flux estimated from Vega spectrum
mag = -2.5*np.log10(totalCheckFlux/filtCorrection)-21.1
print "magnitude of V-filtered spectrum = ", mag


#--------test plots-------------#
ax2 = plt.subplot(1,1,1)
ax2.set_xlabel("Wavelength in $\AA$")
ax2.set_ylabel("Filter Transmission")
plt.plot(filtWvls,filtTrans)
plt.plot(stdWvls,filtTransInterp)
ax2.set_xlim(3000,13000)
plt.show()
#-------------------------------#
#--------test plots-------------#
ax2 = plt.subplot(1,1,1)
ax2.set_xlabel("Wavelength in $\AA$")
ax2.set_ylabel("Filtered Spectrum in W/cm$^{2}$/$\AA$")
plt.plot(stdWvls,stdFlux)
plt.plot(stdWvls,filteredFlux)
ax2.set_xlim(3000,13000)
plt.show()
#-------------------------------#


#multiply filtered spectrum in watts/A by silicon sensor response (A/W)/nm
rname = "221Response.txt"
rdata = np.loadtxt(rname,dtype=float) #load sensor response, units are nm and A/W
respx = np.array(rdata[:,0])*10.0 #convert Wvl to Angstroms
respy = np.array(rdata[:,1])
respInterp = interpolate.griddata(respx,respy,stdWvls,"cubic",fill_value=0)#interpolate response curve to spectrum wvl grid
ampsCurve = respInterp*filteredFlux


#--------test plots-------------#
ax3 = plt.subplot(1,1,1)
ax3.set_xlabel("Wavelength in $\AA$")
ax3.set_ylabel("Optometer Response (A/W)")
plt.plot(respx,respy)
plt.plot(stdWvls,respInterp)
ax3.set_xlim(3000,13000)
plt.show()
#-------------------------------#
#--------test plots-------------#
ax4 = plt.subplot(1,1,1)
ax4.set_xlabel("Wavelength in $\AA$")
ax4.set_ylabel("Converted Spectrum (Amps/$\AA$)")
plt.plot(stdWvls,ampsCurve)
ax3.set_xlim(3000,13000)
plt.show()
#-------------------------------#


#integrate final A/Angs curve to get total Amps expected when viewing Vega through our V filter
totalAmps = integrate.simps(ampsCurve,x=stdWvls)
print "Total Amps expected from Vega V-band = ", totalAmps
print "Measured Amps from optometer given as ", AmpsMeas

#print measured optometer Amps / expected Amps for Vega
print "----------------------------------------"
print "Throughput = measured Amps / expected Amps = \n", AmpsMeas/totalAmps
print "----------------------------------------"





