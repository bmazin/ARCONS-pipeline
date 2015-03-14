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
from scipy.optimize.minpack import curve_fit
from numpy import exp

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
        fraction = 0 #4.0/5.0
        newx = np.arange(int(x[fraction*len(x)]),20000)

        #slopeguess = (np.log(y[-1])-np.log(y[fraction*len(x)]))/(x[-1]-x[fraction*len(x)])
        #print "Guess at exponential slope is %f"%(slopeguess)
        #guess_a, guess_b, guess_c = float(y[fraction*len(x)]), x[fraction*len(x)], slopeguess
        #guess = [guess_a, guess_b, guess_c]

        fitx = x[fraction*len(x)::]
        fity = y[fraction*len(x)::]

        #exp_decay = lambda fx, A, x0, t: A * np.exp((fx-x0) * t)

        #params, cov = curve_fit(exp_decay, fitx, fity, p0=guess, maxfev=2000)
        #A, x0, t= params
        #print "A = %s\nx0 = %s\nt = %s\n"%(A, x0, t)
        #best_fit = lambda fx: A * np.exp((fx-x0)*t)

        #calcx = np.array(newx,dtype=float)
        #newy = best_fit(calcx)

        #normalizing
        norm = fity.max()
        fity/=norm

        guess_a, guess_b = 1/(2*h*c**2/1e-9), 5600 #Constant, Temp
        guess = [guess_a, guess_b]

        blackbody = lambda fx, N, T: N * 2*h*c**2 / (fx)**5 * (exp(h*c/(k*T*(fx))) - 1)**-1 # Planck Law
        #blackbody = lambda fx, N, T: N*2*c*k*T/(fx)**4 #Rayleigh Jeans tail
        #blackbody = lambda fx, N, T: N*2*h*c**2/(fx**5) * exp(-h*c/(k*T*fx)) #Wein Approx

        params, cov = curve_fit(blackbody, fitx*1.0e-8, fity, p0=guess, maxfev=2000)
        N, T= params
        print "N = %s\nT = %s\n"%(N, T)
        best_fit = lambda fx: N * 2*h*c**2 / (fx)**5 * (exp(h*c/(k*T*(fx))) - 1)**-1 #Planck Law
        #best_fit = lambda fx: N*2*c*k*T/(fx)**4 # Rayleigh Jeans Tail
        #best_fit = lambda fx: N*2*h*c**2/(fx**5) * exp(-h*c/(k*T*fx)) #Wein Approx

        calcx = np.array(newx,dtype=float)
        bbfit = best_fit(calcx*1.0E-8)

        calcx = np.array(newx,dtype=float)
        newy = best_fit(calcx*1.0E-8)

        fity*=norm
        newy*=norm

        plt.plot(calcx[3.0*len(fitx)/4.0::],newy[3.0*len(fitx)/4.0::]*1E15,linestyle='--',linewidth=2, color="black",alpha=0.5) #plot fit

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
        plt.plot(re_wl,re_flux*1E15,linestyle="o", marker="o",markersize=6) #plot rebinned spectrum with exp tail
        
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
        f = interpolate.splrep(re_wl,re_flux,w=weight,k=2,s=0)#max(re_flux)**300)
        new_flux = interpolate.splev(new_wl,f,der=0)
        return new_wl, new_flux


c=3.00E10 #cm/s
h=6.626E-27 #erg*s
k=1.3806488E-16 #erg/K

FileName = '/home/srmeeker/scratch/standards/2MASSJ0223_fit_0.npz'
NumFrames = 31

#SETUP PLOTTING
#matplotlib.rcParams.update({'font.size':12, 'font.family': 'sans-serif','sans-serif':['Helvetica']})
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('text',usetex=True)

t = np.load(FileName)

energyBinWidth = 0.1
wvlStart = 3000
wvlStop = 25000
wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth,wvlStart,wvlStop)
nWvlBins = len(wvlBinEdges)-1
binWidths = np.empty(nWvlBins)
for i in xrange(nWvlBins):
    binWidths[i] = wvlBinEdges[i+1]-wvlBinEdges[i]
print binWidths
params = t['params']
wvls = t['wvls']
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]

#print len(wvls)
#print len(binWidths)

curve = 2*np.pi*amps*widths*widths #spectrum of oberved object in counts/s

#print wvls[wvls>4000]
#print wvls[wvls<8000]
#print wvls[(wvls>4000) & (wvls<8000)]
#print curve[wvls>4000]
#print curve[wvls<8000]
totCounts = np.sum(curve[(wvls>4000) & (wvls<25000)])
print "Total counts from 4000 to 8000 Angs = ", totCounts
diam = 5.1055 #5 meter telescope
area = np.pi * ((diam/2.0)**2 -(1.83/2.0)**2) #secondary obstruction diameter 1.83m
print "Total counts per m^2 per second = ", totCounts/area
print "Total counts from Ben's code = 4323"


centers = np.array([3600,4400,5500,6400,7900,12500])
#zeros = [759,1461,999,726,487,194]
#using Ben's method
zeros = np.array([1810,4260,3640,3080,2550,1600])
dlod = np.array([0.15,0.22,0.16,0.23,0.19,0.16])
Jan2Phot = 1.51e7

landVMag = 15.771
#landVCounts = (10.0**(landVMag/(-2.5)))*zeros[2]#*filtCorrectionV # in counts/s/cm^2/Angs
landVCounts = zeros[2]*10.0**(-1*landVMag/2.5)*Jan2Phot*dlod[2]
print "expected V counts from my code= ", landVCounts

print "Throughput = ", ((totCounts/area)/landVCounts)*100, "%"

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(4000,25000)
#ax.set_ylim(0,1E-5)
plt.plot(wvls, curve)
plt.xlabel(ur"Wavelength [\AA]")
plt.ylabel(ur"ARCONS measured Spectrum (Counts/s)")
#plt.savefig("FluxCal_RawSpectrumCounts.eps",format='eps')
plt.show()

#4323 photons per second per m^2 for V=15.771


