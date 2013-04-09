"""
GetFluxCal.py

Created by Seth Meeker on 2012-10-31

Stand-alone FluxCal routine.  Kludged together for preliminary analysis of Lick 2012 Data.

Reads in defocused calibration star pulse-height histogram for each pixel, applies wavelength cal, converts flux from [counts] to [ergs] to match the known spectrum, and normalizes it.
Then reads in the known spectrum of the cal star, in units of [ergs], and already smoothed and normalized,
Then divides the expected spectrum by our measured spectrum for each pixel and produces a scaling relation for each pixel 

All file names and paths are currently hardcoded. Since this is not creating any final products it does not clutter public space on turk. Just change file names in code and run from command line:
> python GetFluxCal.py

"""

import sys
import time
import struct
import os
from os.path import isfile
#import installed libraries
import numpy as np
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from tables import *
from matplotlib.backends.backend_pdf import PdfPages

def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]


herg = 6.6260696E-27 #erg*s
hev = 4.1356675E-15 #eV*s
c = 3.00E18 #Ang/s

roachnumarr = [0,3,4,4,4,4,7,7,1,2,4,4,2,2,2,4,4,4,4,4,4,4,4,4,4]
pixelnumarr = [90,74,7,65,119,175,153,155,149,100,101,149,26,172,0,161,40,52,34,134,179,11,17,163,25]

calpath = "/ScienceData/waveCalSolnFiles/"
fluxpath = "./"
outpath = "/home/srmeeker/FluxCal/SensitivitySpectra/"
if not os.path.exists(outpath):
    os.mkdir(outpath)
    print "making directory for Sensitivity spectra"

junkpath = "/home/srmeeker/scratch/"

#datapath = "/ScienceData/LICK2012/20120914/"
#calfile = "calsol_20120915-050039.h5" #cal for hr7596/58 aquilae
#datafile = "obs_20120915-050547.h5" #data for hr7596/58 aquilae
#skyfile = "obs_20120915-055329.h5" #obs file while searching for css0974 which should be sky
#rawfluxfile = "fhr7596.dat" #real spectrum for hr7596/58 aquilae
#smoothfluxfile = "fhr7596_ergs_smoothed.dat" #smoothed spectrum of hr7596

#datapath = "/ScienceData/LICK2012/20120915/"
#calfile = "calsol_20120916-043616.h5" #cal for hr7596/58 aquilae
#datafile = "obs_20120916-042547.h5" #data for hr7596/58 aquilae
#skyfile =  "obs_20120916-092848.h5" #sky data around Arp 147 observation. Only labeled sky data that night.
#rawfluxfile = "fhr7596.dat" #real spectrum for hr7596/58 aquilae
#smoothfluxfile = "fhr7596_ergs_smoothed.dat" #smoothed spectrum of hr7596
#smoothfluxfile = "fhr7596_bb.dat" #theoretical bb curve for 58 Aqu

datapath = "/ScienceData/LICK2012/20120915/"
calfile = "calsol_20120916-043616.h5" #cal for hd176670
datafile = "obs_20120916-040910.h5" #data for hd176670
skyfile =  "obs_20120916-092848.h5" #sky data around Arp 147 observation. Only labeled sky data that night.
rawfluxfile = "HD176670_bb.dat" #real spectrum for hd176670 aquilae
smoothfluxfile = "HD176670_bb.dat" #smoothed spectrum of hd176670

#datapath = "/ScienceData/LICK2012/20120916/"
#calfile = "calsol_20120917-071245.h5" #cal file for BD+25 4655
#datafile = "obs_20120917-072009.h5" #data file for BD+25 4655
#rawfluxfile = "fbd25d4655_ergs.dat" #real spectrum for BD+25 4655
#smoothfluxfile = "fbd25d4655_smoothed.dat" #smoothed spectrum for BD+25

junkfile = "junkplots.ps"
outfile = "FluxCal.h5"
outtxt = "Sensitivity.txt"

#leave this in for debugging.  Eventually take it out and only supply the smoothed calibration spectrum
#read in expected spectrum
fname = fluxpath + rawfluxfile
fdata = np.loadtxt(fname,dtype=float)
rawwl = np.array(fdata[::,0])
rawflux = fdata[::,1]/(1E16)
rawflux = np.array(rawflux)
#convert rawflux in ergs to counts
E2n = rawwl/(herg*c)
rawflux = rawflux*E2n
#normalize
norm = float(rawflux.max())
#rawflux = rawflux/norm

#THIS IS NOW DONE IN ANOTHER ROUTINE AND SMOOTHED DATA IS PASSED IN TO HERE
#smooth real spectrum to get rid of absorption features
#wlen=len(inwl)/20 #window length
#flux = smooth(rawflux,window_len=wlen,window='flat')
#normalize
#influx = flux/flux.max()

fname = fluxpath + smoothfluxfile
fdata = np.loadtxt(fname,dtype=float)
inwl = np.array(fdata[::,0])
influx = fdata[::,1]
norm = influx.max()
influx = influx/norm

#read in energy scaling file
try:
    calh5 = openFile(str(calpath+calfile),'r')
except:
    print "unable to open wavelength cal file"
    sys.exit()
caltable = calh5.root.wavecal.calsoln.read()
calh5.close()

#read in our calibration star data
try:
    h5file = openFile(str(datapath+datafile), 'r')
except:
    print "unable to open calibration star observation file"
    sys.exit()

#read in sky data
try:
    skyfile = openFile(str(datapath+skyfile),'r')
except:
    print "unable to open sky file for flux cal"
    sys.exit()

#get sky exposure time for subtraction later and sky bmap
skybmap = skyfile.root.beammap.beamimage.read()
htable = skyfile.root.header.header.read()
tsky = int(htable["exptime"])

#get data exposure time for subtraction later
bmap = h5file.root.beammap.beamimage.read()
htable = h5file.root.header.header.read()
tobs = int(htable["exptime"])
nxpix = np.shape(bmap)[1]
nypix = np.shape(bmap)[0]
print "nxpix= ",nxpix
print "nypix= ",nypix

#for n in xrange(len(roachnumarr)):
#    roachnum = roachnumarr[n]
#    pixelnum = pixelnumarr[n]
#    print roachnum, pixelnum

for l in xrange(len(caltable)):
    roachnum = caltable[l][4]
    pixelnum = caltable[l][1]
    print "\n-----------------------------\n"
    print l
    print "cal pixel /r%s/p%s"%(roachnum,pixelnum)
    x = caltable[l][0]
    y = caltable[l][2]
    xoff = caltable[l][3][0]
    yoff = caltable[l][3][1]
    amp = caltable[l][3][2]
    flag = caltable[l][7]
    #print xoff, amp, yoff
    print "(%i,%i)"%(x,y)
    print "flag = ", flag
    if flag != 0:
        print "BAD CALIBRATION"
        continue
    
    #read in sky data first
    photons= np.concatenate(skyfile.root._f_getChild(skybmap[y][x])[:])
    peakheights= np.right_shift(photons,32)%4096
    parabheights= np.right_shift(photons,44)%4096
    baseline= np.right_shift(photons,20)%4096
    obsheights = np.array(parabheights,dtype=float)-np.array(baseline,dtype=float)
    #apply scaling polynomial to convert pulse heights to wavelength in Angstroms
    ev = amp*(obsheights-xoff)**2 + yoff
    for v in xrange(len(ev)):
        if ev[v] <= 0.5:
            ev[v] = 0.5 #check for 0 energy photons which will fail division later
    lambdas = (hev*c)/ev
    reallambdas = np.where(lambdas>=1000) #cut out photons less than 1000 Angstroms
    lambdas = lambdas[reallambdas]
    reallambdas = np.where(lambdas<=20000) #cut out photons greater than 20,000 Angs
    lambdas = lambdas[reallambdas]
    skylambdas = lambdas

    #then read in obs data (this should really be its own function that gets called)
    print "bmap pixel", (str(bmap[y][x]))
    photons= np.concatenate(h5file.root._f_getChild(bmap[y][x])[:])
    peakheights= np.right_shift(photons,32)%4096
    parabheights= np.right_shift(photons,44)%4096
    baseline= np.right_shift(photons,20)%4096
    obsheights = np.array(parabheights,dtype=float)-np.array(baseline,dtype=float)
    #apply scaling polynomial to convert pulse heights to wavelength in Angstroms
    ev = amp*(obsheights-xoff)**2 + yoff
    for v in xrange(len(ev)):
        if ev[v] <= 0.5:
            ev[v] = 0.5 #check for 0 energy photons which will fail division later
    lambdas = (hev*c)/ev
    reallambdas = np.where(lambdas>=1000) #cut out photons less than 1000 Angstroms
    lambdas = lambdas[reallambdas]
    reallambdas = np.where(lambdas<=20000) #cut out photons greater than 20,000 Angs
    lambdas = lambdas[reallambdas]

    if len(lambdas) == 0 or len(skylambdas) == 0:
        print "BAD CAL/WRONG FLAG: 1ST WL-CUT REMOVED ALL PHOTONS"
        print "skipping...."
        continue

    #make histogram of obsheights
    nobsbins = obsheights.max() - obsheights.min()
    obshist,obsbins = np.histogram(obsheights,bins=nobsbins,range=(obsheights.min(),obsheights.max()))
    obshist = obshist/float(obshist.max())

    #make histogram of lambdas, now a spectrum in flux units of [counts]
    nlambins = abs(lambdas.max() - lambdas.min())
    try:
        rawlamhist,rawlambins = np.histogram(lambdas,bins=nlambins,range=(lambdas.min(),lambdas.max()))
        indices = np.where(rawlamhist != 0)
        rawlamhist = rawlamhist[indices]
        rawlambins = rawlambins[indices]
    except ValueError:
        print "BAD CAL/WRONG FLAG: ONLY 1 HISTOGRAM BIN"
        print "skipping...."
        continue

    #make histogram of sky photons, binned to same histogram as real data
    try:
        rawskyhist,rawskybins = np.histogram(skylambdas,bins=nlambins,range=(lambdas.min(),lambdas.max()))
        indices = np.where(rawskyhist != 0)
        rawskyhist = rawskyhist[indices]
        rawskybins = rawskybins[indices]
        #put sky histogram onto same grid as data
        rawskyhist = np.interp(rawlambins,rawskybins,rawskyhist)
    except ValueError:
        print "BAD CAL/WRONG FLAG: ONLY 1 SKY HISTOGRAM BIN"
        print "skipping...."
        continue

    #convert from histogram of flux in counts to flux in ergs
    n2erg = herg*c/rawlambins
    rawlamhist = rawlamhist*n2erg
    rawskyhist = rawskyhist*n2erg

    #convert both data and sky fluxes into [ergs/s] and subtract scaled skyhistogram from raw data histogram
    rawlamhist = rawlamhist/float(tobs)
    rawskyhist = rawskyhist/float(tsky)
    skysubtracted = rawlamhist-rawskyhist #now in units of [ergs/s/Angs]

    #use this when skysubtraction is turned off
    #skysubtracted = rawlamhist

    #smooth raw spectrum
    wlen=len(rawlamhist)/50 #window length
    rawsmooth = smooth(rawlamhist,window_len=wlen,window='flat')

    #smooth sky subtracted spectrum
    wlen=len(skysubtracted)/50 #window length
    lamhist = smooth(skysubtracted,window_len=wlen,window='flat')

    #normalize
    lamhist = lamhist/lamhist.max()
    rawlamhist = rawlamhist/rawlamhist.max()
    rawskyhist = rawskyhist/rawlamhist.max() #normalize skyhist to obs hist to subtract the right fraction of sky counts

    #make test plots
    pp = PdfPages(outpath+'r%i_p%i_plots.pdf'%(roachnum, pixelnum))
    matplotlib.rcParams['font.size']=6

    plt.figure()
    #ax1 = plt.subplot(321)
    #ax1.set_title('obshist')
    #plt.plot(obsbins[:-1],obshist)
    #plt.show()
    ax2 = plt.subplot(332)
    ax2.set_title('our flux in ergs/s')
    plt.plot(rawlambins,rawlamhist)
    #plt.show()
    ax3 = plt.subplot(331)
    ax3.set_title('real spectrum')
    plt.plot(rawwl,rawflux)
    #plt.show()
    ax4 = plt.subplot(333)
    ax4.set_title('raw sky spectrum')
    plt.plot(rawlambins,rawskyhist)

    #cut spectra to same wavelength range

    lmin = np.max([rawlambins.min(),inwl.min()])
    lmax = np.min([rawlambins.max(),inwl.max()])

    try:
        indmin = np.where(inwl >= lmin)[0][0]
        indmax = np.where(inwl <= lmax)[0][-1]
        wl = inwl[indmin:indmax]
        flux = influx[indmin:indmax]

        indmin = np.where(rawlambins >= lmin)[0][0]
        indmax = np.where(rawlambins <= lmax)[0][-1]
        lambins = rawlambins[indmin:indmax]
        rawsmooth = rawsmooth[indmin:indmax]
        lamhist = lamhist[indmin:indmax]
    except IndexError:
        print "BAD CAL/WRONG FLAG: DATA WL RANGE OUTSIDE FLUX CAL WL RANGE"
        print "skipping...."
        pp.savefig()
        pp.close()
        continue

    #more test plots
    ax11 = plt.subplot(334)
    ax11.set_title('raw counts minus sky')
    plt.plot(rawlambins,skysubtracted)
    ax12 = plt.subplot(335)
    ax12.set_title('cut and smoothed counts')
    plt.plot(lambins,rawsmooth)
    ax5 = plt.subplot(337)
    ax5.set_title('cut and smoothed counts minus sky')
    plt.plot(lambins,lamhist)
    #plt.show()
    ax6 = plt.subplot(336)
    ax6.set_title('wl cut real spectrum')
    plt.plot(wl,flux)
    #plt.show()

    if len(lambins) == 0:
        print "BAD CAL/WRONG FLAG: 2ND WL-CUT REMOVED ALL PHOTONS"
        print "skipping...."
        pp.savefig()
        pp.close()
        continue

    #interpolate measured spectrum to same wl grid as real spectrum
    newgrid = xrange(int(lmax-lmin)) + lmin
    newflux = np.interp(newgrid,wl,flux)
    newhist = np.interp(newgrid,lambins,lamhist)

    #divide real spectrum by measured spectrum to get "sensitivity spectrum"
    sensitivity = newflux/newhist

    #set any negative factors to 1
    sensitivity[sensitivity <= 0] = 1.0
    sensitivity[sensitivity == np.inf] = 1.0

    #more test plots
    #ax7 = plt.subplot(325)
    #ax7.set_title('interp\'d lambdas')
    #plt.plot(newgrid,newhist)
    #plt.show()
    #ax8 = plt.subplot(326)
    #ax8.set_title('interp\'d spectrum')
    #plt.plot(newgrid,newflux)
    #plt.show()
    ax9 = plt.subplot(338)
    ax9.set_title('sensitivity spectrum')
    plt.plot(newgrid,sensitivity)
    #plt.show()

    #check flux calibration on original data
    newrawhist = np.interp(newgrid,rawlambins[indmin:indmax],rawlamhist[indmin:indmax])
    fluxcalibrated = sensitivity*newrawhist

    ax10 = plt.subplot(339)
    ax10.set_title('flux calibrated spectrum')
    plt.plot(newgrid, fluxcalibrated)
    #plt.show()

    outarr = np.empty((len(newgrid),2),dtype=float)
    outarr[:,0]=newgrid
    outarr[:,1]=sensitivity
    #save sensitivity spectrum to file
    fname = str(outpath+'r%i_p%i_sensitivity.txt'%(roachnum,pixelnum))
    np.savetxt(fname, outarr)


    #plt.legend()
    pp.savefig()
    pp.close()

h5file.close()
skyfile.close()

print "\n-------------------\n"
print "FINISHED FLUX CAL WITH ", datafile
print "\n--------------------\n"


