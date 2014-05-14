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
Written by Seth Meeker, 11/20/13
Given a B, V, or R magnitude,
an ARCONS QE file, 
and a measured telescope throughput in that band:
spits out how many counts to expect on the array

Calling sequence:
> python expectedCounts.py [filter band] [magnitude] [telescope efficiency]
Example:
> python expectedCounts.py V 14.9 0.239

"""

if len(sys.argv)<3:
    print "ERROR: NOT ENOUGH ARGUMENTS IN CALLING SEQUENCE"
    print "Provide desired band, star magnitude, and measured telescope efficiency in that band."
    print "Example: >python expectedCounts.py V 14.9 0.239"
    print "Available filters: B, V, R"
    print "If no efficiency is provided, the hardcoded efficiency for the given band will be used\n----------------------\n"
    sys.exit()

band = str(sys.argv[1])
mag = float(sys.argv[2])

QEFileName = "ARCONS_QE_2013_fake.txt"

if band == ('b' or 'B'):
    fwhm = 940
    lambda_c = 4400
    dlam_over_lam = 0.22
    flux_0 = 4260 #janskys
    scaling = 0.9871/0.7300
    eff = 0.20 #hardcoded B band telescope throughput

if band == ('v' or 'V'):
    fwhm = 880
    lambda_c = 5500
    dlam_over_lam = 0.16
    flux_0 = 3640 #janskys
    scaling = 0.9799/0.8761
    eff = 0.20 #hardcoded v band telescope throughput

if band == ('r' or 'R'):
    fwhm = 1380
    lambda_c = 6400
    dlam_over_lam = 0.23
    flux_0 = 3080 #janskys
    scaling = 0.9838/0.8319
    eff = 0.20 #hardcoded R band telescope throughput

if len(sys.argv)>3:
    eff = float(sys.argv[3]) #if a telescope throughput is provided in the arguments, use that

##!!!!!!!! NEED TO CHANGE FLUX_0 FOR EACH BAND SINCE OUR BVR FILTERS ARE NOT TYPICAL COLORED GLASS !!!!!!!!!
# multiply flux_0 by a scaling factor assuming typical maximum J-C filter transmissions for the quoted values
# and using the maximum measured values for our own J-C filters
flux_0 *= scaling

#convert mag to counts using conversions given above. Don't forget to multiply by telescope collecting area.
diam = 5.1055 #5 meter telescope
area = np.pi * ((diam/2.0)**2 -(1.83/2.0)**2) #secondary obstruction diameter 1.83m

photonsPerSecPerMeterSquared = flux_0 * 10.0**(-mag/2.5) * dlam_over_lam *1.51E7 #1.51e7/(dlam/lam) photons/s/m2 = 1 Jy
p_per_s = photonsPerSecPerMeterSquared * area #photons/s collected by telescope
print "\n#-------------------------------------#\nGiven %s %s magnitude\n%s photons per second collected by telescope."%(mag, band,int(p_per_s))

#multiply by telescope efficiency to get counts per second at ARCONS
p_per_s *= eff
print "\nGiven %s %% efficiency from the telescope\n%s photons per second make it to ARCONS."%(eff*100,int(p_per_s))

#find ARCONS QE. Load QE file and find %QE average of requested band.
QEfile = os.environ['ARCONS_PIPELINE_PATH']+'/util/data/'+QEFileName
fdata = np.loadtxt(QEfile,dtype=float)
wvls = np.array(fdata[:,0])*10.0 #convert from nm to Angstroms
QEcurve = np.array(fdata[:,1])
start = lambda_c - fwhm/2.0
stop = lambda_c + fwhm/2.0
qe = np.mean(QEcurve[(wvls>start) & (wvls<stop)])

#multiply counts calculated by ARCONS QE and telescope efficiency to get final expected counts in that band
p_per_s *= qe
print "\nWith %s %% QE from ARCONS\n%s photons per second should be measured by ARCONS.\n#-------------------------------------#\n"%(qe*100,int(p_per_s))










