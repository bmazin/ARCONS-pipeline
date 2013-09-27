from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from util.rebin import rebin
import matplotlib
from scipy import interpolate
from scipy import integrate
from scipy.optimize.minpack import curve_fit
from math import exp


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
    print (filty)

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


c=3.00E18 #Angs/s
h=6.626E-27 #erg*s

#Vfile = '/home/srmeeker/scratch/standards/Landolt9542_V_fit.npz'
#Rfile = '/home/srmeeker/scratch/standards/Landolt9542_R_fit.npz'
Vfile = '/home/srmeeker/ARCONS-pipeline/examples/Pal2012_landoltPhot/Landolt9542_wOFF_hpON/Landolt9542_V_raw_0.npz'
Rfile = '/home/srmeeker/ARCONS-pipeline/examples/Pal2012_landoltPhot/Landolt9542_wOFF_hpON/Landolt9542_R_raw_0.npz'
NumFrames = 31

vspec = np.zeros((NumFrames))
rspec = np.zeros((NumFrames))
vf = np.load(Vfile)
rf = np.load(Rfile)

energyBinWidth = 0.1
wvlStart = 3000
wvlStop = 13000
wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth,wvlStart,wvlStop)
nWvlBins = len(wvlBinEdges)-1
binWidths = np.empty(nWvlBins)
for i in xrange(nWvlBins):
    binWidths[i] = wvlBinEdges[i+1]-wvlBinEdges[i]
#print binWidths

vstack = vf['stack']
rstack = rf['stack']
wvls = vf['wvls']
print wvls

#print len(wvls)
#print len(binWidths)

for i in xrange(NumFrames):
    vframe = vstack[i]
    rframe = rstack[i]
    vframe[np.isnan(vframe)]=0
    rframe[np.isnan(rframe)]=0
    vspec[i] = np.sum(vframe)
    rspec[i] = np.sum(rframe)

print vspec
print rspec

#plot v-band filter for reference
fname = 'Vfilter.txt'
#fname = 'Subaru_V-Band_Transmission.txt'
tx,ty, junk, junk = loadFilter(fname)
ty = ty/ max(ty) * .8623 #normalize filter to our Johnson V filter which has max transmission of 86.23%

#plot v-band filter for reference
fname = 'Rfilter.txt'
#fname = 'Subaru_V-Band_Transmission.txt'
txR,tyR, junk, junk = loadFilter(fname)
txR = txR*10.0
tyR = tyR/max(tyR)*0.8218 #normalize filter to our Johnson R filter which has max transmission of 82.18%

#plot v-band filter for reference
fname = 'Subaru_R-Band_Transmission.txt'
txSR, tySR, junk, junk = loadFilter(fname)

plt.plot(tx,ty,'green')
plt.plot(txR, tyR,'red')
plt.legend(["Johnson V Filter","Johnson R Filter"])
plt.title("ARCONS Filters")
print max(ty)
plt.show()

plt.plot(wvls,vspec,color='green')
plt.plot(wvls,rspec,color='red')
plt.plot(tx,ty/max(ty)*max(vspec[:-len(vspec)/3]),color = 'black',alpha =0.4)
plt.plot(txR,tyR/max(tyR)*max(rspec[:-len(rspec)/4]),color = 'black',alpha =0.7)
#plt.plot(txSR,tySR/max(tySR)*max(rspec[:-len(rspec)/4]),color = 'black',alpha =1.0)
plt.show()
    

