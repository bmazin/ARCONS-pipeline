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
        #plot this to debug std spectrum fitting
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
        weight = [1,1,1,1,1,1,1]
        #print len(weight)
        #weight = re_flux/min(re_flux)
        #weight = 1.0/weight
        #weight = weight/max(weight)
        print weight
        f = interpolate.splrep(re_wl,re_flux,w=weight,k=3,s=max(re_flux)**200)
        new_flux = interpolate.splev(new_wl,f,der=0)
        return new_wl, new_flux


c=3.00E18 #Angs/s
h=6.626E-27 #erg*s

NumFrames = 31
nFiles = 13
curves = np.zeros((nFiles,NumFrames),dtype=float)

for k in xrange(nFiles):
    FileName = '/home/srmeeker/scratch/standards/sdssj0926_fit_%s.npz'%(k)
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
print "Median Spectrum = "
print curve

#SETUP PLOTTING
fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlim(3700,13000)
#plt.ylim(0,0.001)

#LOAD KNOWN SPECTRUM OF STANDARD
objectName = "sdss j0926"
'''
#import the known spectrum of the calibrator and rebin to the histogram parameters given
#must be imported into array with dtype float so division later does not have error
std = MKIDStd.MKIDStd()
a = std.load(objectName)
x = a[:,0]
y = np.array(a[:,1]) #std object spectrum in ergs/s/Angs/cm^2
#print y
y*=(1E-17) #SDSS J0926 saved in units of 10^-17 ergs/s/Angs/cm^2
print y
#convert from ergs/s to counts/s
'''

#Manually load j0926 standard file. MKIDStd units may be wrong.
scale = (1.0E17) #for .txt files strange unit scaling
fname = 'sdssj0926.txt'
fdata = np.loadtxt(fname,dtype=float)
x = np.array(fdata[:,0])
y = np.array(fdata[:,1])/scale #flux in ergs/s/cm^2/Angs
print y
y = y/(h*c/x) #flux in counts/s/cm^2/Angs
print y

print "Loaded standard spectrum of %s"%(objectName)

newwl, newflux = cleanSpectrum(x,y,objectName,wvlBinEdges)
print newwl
print newflux
#plot these to debug std spectrum fitting
#plt.plot(newwl,newflux)
#plt.plot(x,y)
#plt.show()

newa = rebin(newwl,newflux,wvlBinEdges)
x = newa[:,0]
y = newa[:,1]

#print x
#print y

ax.plot(wvls,curve)
ax.plot(x,y)
ax.set_title('Absolute Spectrum of  '+FileName.split('/')[-1].split('_')[0])
ax.set_yscale('log')

plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wvls,curve/y)
ax.set_ylim(0,0.06)
#ax.set_xlim(3750,12500)
plt.show()

np.savez('%s_throughput.npz'%(objectName),throughput=curve/y,wvls=wvls)

