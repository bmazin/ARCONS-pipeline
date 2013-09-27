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

def cleanSpectrum(x,y,objectName, wvlBinEdges):
        #locations and widths of absorption features in Angstroms
        #features = [3890,3970,4099,4340,4860,6564,6883,7619]
        #widths = [50,50,50,50,50,50,50,50]
        #for i in xrange(len(features)):
        #    #check for absorption feature in std spectrum
        #    ind = np.where((x<(features[i]+15)) & (x>(features[i]-15)))[0]
        #    if len(ind)!=0:
        #        ind = ind[len(ind)/2]
        #    #if feature is found (flux is higher on both sides of the specified wavelength where the feature should be)
        #    if y[ind]<y[ind+1] and y[ind]<y[ind-1]:
        #        #cut out width[i] around feature[i]
        #        inds = np.where((x >= features[i]+widths[i]) | (x <= features[i]-widths[i]))
        #        x = x[inds]
        #        y = y[inds]

        #fit a tail to the end of the spectrum to interpolate out to desired wavelength in angstroms
        fraction = 2.1/3
        newx = np.arange(int(x[fraction*len(x)]),20000)

        slopeguess = (np.log(y[-1])-np.log(y[fraction*len(x)]))/(x[-1]-x[fraction*len(x)])
        print "Guess at exponential slope is %f"%(slopeguess)
        guess_a, guess_b, guess_c = float(y[fraction*len(x)]), x[fraction*len(x)], slopeguess
        guess = [guess_a, guess_b, guess_c]

        fitx = x[fraction*len(x):]
        fity = y[fraction*len(x):]

        exp_decay = lambda fx, A, x0, t: A * np.exp((fx-x0) * t)

        params, cov = curve_fit(exp_decay, fitx, fity, p0=guess, maxfev=2000)
        A, x0, t= params
        print "A = %s\nx0 = %s\nt = %s\n"%(A, x0, t)
        best_fit = lambda fx: A * np.exp((fx-x0)*t)

        calcx = np.array(newx,dtype=float)
        newy = best_fit(calcx)

        #func = interpolate.splrep(x[fration*len(x):],y[fraction*len(x):],s=smooth)
        #newx = np.arange(int(x[fraction*len(x)]),self.wvlBinEdges[-1])
        #newy = interpolate.splev(newx,func)

        wl = np.concatenate((x,newx[newx>max(x)]))
        flux = np.concatenate((y,newy[newx>max(x)]))

        #new method, rebin data to grid of wavelengths generated from a grid of evenly spaced energy bins
        #R=7.0 at 4500
        #R=E/dE -> dE = R/E
        dE = 0.3936 #eV
        start = 1000 #Angs
        stop = 25000 #Angs
        enBins = ObsFile.makeWvlBins(dE,start,stop)
        rebinned = rebin(wl,flux,enBins)
        re_wl = rebinned[:,0]
        re_flux = rebinned[:,1]
        #plt.plot(re_wl,re_flux,color='r')
        
        re_wl = re_wl[np.isnan(re_flux)==False]
        re_flux = re_flux[np.isnan(re_flux)==False]

        start1 = wvlBinEdges[0]
        stop1 = wvlBinEdges[-1]
        #regrid downsampled data 

        new_wl = np.arange(start1,stop1)

        #print re_wl
        #print re_flux
        #print new_wl

        #weight=1.0/(re_flux)**(2/1.00)
        print len(re_flux)
        weight = np.ones(len(re_flux))
        #decrease weights near peak
        ind = np.where(re_flux == max(re_flux))[0]
        weight[ind] = 0.3
        for p in [1,2,3]:
            if p==1:
                wt = 0.3
            elif p==2:
                wt = 0.6
            elif p==3:
                wt = 0.7
            try:
                weight[ind+p] = wt
            except IndexError:
                 pass
            try:
                 if ind-p >= 0:
                     weight[ind-p] = wt
            except IndexError:
                pass
        #change weights to set how tightly fit must match data points
        #weight[-4:] = 1.0
        weight = [0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7]
        #print len(weight)
        #weight = re_flux/min(re_flux)
        #weight = 1.0/weight
        #weight = weight/max(weight)
        print weight
        f = interpolate.splrep(re_wl,re_flux,w=weight,k=3,s=max(re_flux)**200)
        new_flux = interpolate.splev(new_wl,f,der=0)
        return new_wl, new_flux

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
"""
#LOAD ALL SDSSj0926/CRAB FILES

NumFrames = 31
nFiles = 13
curves = np.zeros((nFiles,NumFrames),dtype=float)

for k in xrange(nFiles):
    FileName = '/home/srmeeker/scratch/standards/sdssj0926/sdssj0926_fit_%s.npz'%(k)
    print FileName
    t = np.load(FileName)

    energyBinWidth = 0.1
    wvlStart = 3000
    wvlStop = 13000
    wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth,wvlStart,wvlStop)
    nWvlBins = len(wvlBinEdges)-1
    binWidths = np.empty(nWvlBins)
    for i in xrange(nWvlBins):
        binWidths[i] = wvlBinEdges[i+1]-wvlBinEdges[i]
    #print binWidths
    params = t['params']
    wvls = t['wvls']
    amps = params[:,1]
    widths = params[:,4]
    xpos = params[:,2]
    ypos = params[:,3]

    #print len(wvls)
    #print len(binWidths)

    curve = 2*np.pi*amps*(widths**2) #spectrum of oberved object in counts/s
    curve /= binWidths #spectrum is now in counts/s/Angs

    diam = 500 #5 meter telescope
    area = np.pi * ((diam/2.0)**2 -(diam/4.0)**2)
    curve/= area #spectrum is now in counts/s/Angs/cm^2
    print k
    print curve
    curves[k] = curve
print curves
finalCurve = np.zeros((NumFrames),dtype=float)
for j in xrange(NumFrames):
    finalCurve[j] = np.median(curves[:,j])

curve = finalCurve

"""
#set up magnitude parameters
# bands are [U,B,V,R,I,J]
#mags = [16.68,17.22,16.64,16.14,15.61,14.72] #from Sandberg and Sollerman, 2009, Crab+Knot mags
errors = [0.03,0.02,0.02,0.01,0.01,0.03]

centers = np.array([3600,4400,5500,6400,7900,12500])
#zeros = [759,1461,999,726,487,194]

#using Ben's method
zeros = np.array([1810,4260,3640,3080,2550,1600])
dlod = np.array([0.15,0.22,0.16,0.23,0.19,0.16])
Jan2Phot = 1.51e7

#energyBinWidth = 0.3936
energyBinWidth = 0.1
wvlStart = 3000
wvlStop = 13000
wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth,wvlStart,wvlStop)
nWvlBins = len(wvlBinEdges)-1
binWidths = np.empty(nWvlBins)
for i in xrange(nWvlBins):
    binWidths[i] = wvlBinEdges[i+1]-wvlBinEdges[i]
#print binWidths

#LOAD A SINGLE FILE
#FileName = '/home/srmeeker/scratch/standards/crab_fit_11.npz'
#FileName = '/home/srmeeker/scratch/standards/sdssj0926/sdssj0926_fit_6.npz'
#FileName = '/home/srmeeker/scratch/standards/corot18/corot18_fit.npz'
#NumFrames = 31

#GET TOTAL VBAND COUNTS
FileName = '/home/srmeeker/scratch/standards/Landolt9542_V_fit.npz'
t = np.load(FileName)
#LOAD GAUSSIAN FIT PARAMETERS
params = t['params']
wvls = t['wvls']
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
#INTEGRATE UNDER GAUSSIAN TO GET TOTAL COUNTS IN PSF
totalVCounts = 2*np.pi*amps*(widths**2) #spectrum of observed object in counts/s
print "totalVCounts for object = ",totalVCounts



#GET TOTAL RBAND COUNTS
FileName = '/home/srmeeker/scratch/standards/Landolt9542_R_fit.npz'
t = np.load(FileName)
#LOAD GAUSSIAN FIT PARAMETERS
Rparams = t['params']
wvls = t['wvls']
amps = Rparams[:,1]
widths = Rparams[:,4]
xpos = Rparams[:,2]
ypos = Rparams[:,3]
#INTEGRATE UNDER GAUSSIAN TO GET TOTAL COUNTS IN PSF
totalRCounts = 2*np.pi*amps*(widths**2) #spectrum of observed object in counts/s
print "totalRCounts for object = ",totalRCounts


#telescope parameters
diam = 500 #5 meter telescope
area = np.pi * ((diam/2.0)**2 -(100)**2) #Palomar secondary is ~1m radius from Serabyn 2007

#convert area into m^2
area*= 1e-4

totalVCounts /= area # counts/s/m^2
totalRCounts /= area

print "VCounts/s/m^2", totalVCounts
print "RCounts/s/m^2", totalRCounts

'''
Vfname = 'Vfilter.txt'
Vx,Vy,filtWidthV,filtCorrectionV = loadFilter(Vfname)
Rfname = 'Rfilter.txt'
Rx,Ry,filtWidthR,filtCorrectionR = loadFilter(Rfname)

totalVCounts = totalVCounts#/(filtWidthV) #divide by width of filter to get counts/s/cm^2/Angstrom
totalRCounts = totalRCounts#/(filtWidthR)
'''

#landolt9542 parameters
landVMag = 15.61
#landVCounts = (10.0**(landVMag/(-2.5)))*zeros[2]#*filtCorrectionV # in counts/s/cm^2/Angs
landVCounts = zeros[2]*10.0**(-1*landVMag/2.5)*Jan2Phot*dlod[2]
print "landolt9542 expected V counts = ", landVCounts
print "landolt9542 measured V counts = ", totalVCounts
landVTP = totalVCounts/landVCounts #Vband

#V-R = -0.119
#R = V+0.119

landRMag = landVMag+0.119
#landRCounts = (10.0**(landRMag/(-2.5)))*zeros[3]#*filtCorrectionR # in counts/s/cm^2/Angs
landRCounts = zeros[3]*10.0**(-1*landRMag/2.5)*Jan2Phot*dlod[3]
print "landolt9542 expected R counts = ", landRCounts
print "landolt9542 measured R counts = ", totalRCounts
landRTP = totalRCounts/landRCounts #Vband

print "Landolt V-band TP = ",landVTP
print "Landolt R-band TP = ",landRTP

#SETUP PLOTTING
fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlim(4000,12500)
plt.ylim(5e-3,6e-2)

fname = 'throughput.txt'
fdata = np.loadtxt(fname,dtype=float)
xtp = np.array(fdata[:,0])
tp = np.array(fdata[:,1])

fname = 'throughput_high.txt'
fdata = np.loadtxt(fname,dtype=float)
xtph = np.array(fdata[:,0])
tph = np.array(fdata[:,1])

plt.plot(xtp,tp,color="black")
plt.plot(xtph,tph,color="grey")
plt.plot(centers[2],landVTP,'o',color="green")
plt.plot(centers[3],landRTP,'o',color="red")

#ax.set_yscale('log')

plt.show()

