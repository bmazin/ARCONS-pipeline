from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from headers import pipelineFlags
from util.utils import gaussianConvolution
import matplotlib
from mpltools	import style
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
        fraction = 4.0/5.0
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

        plt.plot(calcx,newy) #plot exponential tail

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
        plt.plot(re_wl,re_flux,) #plot rebinned spectrum with exp tail
        
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

##=======================  Define some Constants     ============================
c=3.00E10 #cm/s
h=6.626E-27 #erg*s
k=1.3806488E-16 #erg/K
heV = 4.13566751E-15


##==========================   Import Std Spectrum   ============================
objectName = "G158-100"
#import the known spectrum of the calibrator and rebin to the histogram parameters given
#must be imported into array with dtype float so division later does not have error
std = MKIDStd.MKIDStd()
a = std.load(objectName)
a = std.countsToErgs(a)
x = a[:,0]
y = np.array(a[:,1]) #std object spectrum in ergs/s/Angs/cm^2
plt.plot(x,y)
plt.show()


##=============== BB Fit to extend spectrum to 11000 Angstroms ==================
fraction = 1.0/5.0 #how much of original spectrum do you want to EXCLUDE from the fit?
newx = np.arange(int(x[fraction*len(x)]),20000)
fitx = np.array(x[int(fraction*len(x))::])
fity = np.array(y[int(fraction*len(x))::])
norm = fity.max()
fity/=norm

guess_a, guess_b = 1/(2*h*c**2/1e-9), 5600 #Constant, Temp
guess = [guess_a, guess_b]

blackbody = lambda fx, N, T: N * 2*h*c**2 / (fx)**5 * (np.exp(h*c/(k*T*(fx))) - 1)**-1 # Planck Law
params, cov = curve_fit(blackbody, fitx*1.0e-8, fity, p0=guess, maxfev=2000)
N, T= params

print "N = %s\nT = %s\n"%(N, T)
best_fit = lambda fx: N * 2*h*c**2 / (fx)**5 * (np.exp(h*c/(k*T*(fx))) - 1)**-1 #Planck Law
calcx = np.array(newx,dtype=float)
bbfit = best_fit(calcx*1.0E-8)

calcx = np.array(newx,dtype=float)
newy = best_fit(calcx*1.0E-8)
fity*=norm
newy*=norm

plt.plot(x, y)
plt.plot(calcx,newy,linestyle='--',linewidth=2, color="black",alpha=0.5) #plot fit
plt.show()

wl = np.concatenate((x,newx[newx>max(x)]))
flux = np.concatenate((y,newy[newx>max(x)]))
x=wl
y=flux

plt.plot(x, y)
plt.show()

'''
##======  multiply spectrum by ARCONS filter transmission to properly simulate edge cutoffs ==
transmission = np.ones(len(x))
transmission[x<3500]=0
transmission[x>11500]=0
y*=transmission
plt.plot(x,transmission)
plt.show()
plt.plot(x,y)
plt.show()
'''

newX, newY = gaussianConvolution(x,y,xEnMin=0.005,xEnMax=6.0,xdE=0.001,fluxUnits = "lambda",r=8,plots=True)



'''
##================  Convert to F_nu and put x-axis in frequency  ===================
xEn = heV*(c*1E8)/x
xNu = xEn/heV
yNu = y * x**2 * 3.34E4 #convert Flambda to Fnu(Jy)
plt.plot(xNu,yNu)
plt.show()



##============  regrid to a constant energy spacing for convolution  ===============
xNuGrid = np.arange(0.005,6.0,0.001)/heV #make new x-axis gridding in constant freq bins
yNuGrid = interpolate.griddata(xNu,yNu,xNuGrid, 'linear', fill_value=0)
plt.plot(xNuGrid,yNuGrid)
#plt.plot(xNuGrid*heV, yNuGrid/heV)
#plt.show()



#integrate to get total flux so we can make sure flux is conserved after convolution ==



##======  define gaussian for convolution, on same gridding as spectral data  ======
#WARNING: right now flux is NOT conserved
amp = 1.0
offset = 0
dE = 0.4
sig = dE/heV/2.355 #define sigma as FWHM (.4eV) converted to frequency
gaussX = np.arange(-2,2,0.001)/heV
gaussY = amp * np.exp(-1.0*(gaussX-offset)**2/(2.0*(sig**2)))
plt.plot(gaussX, gaussY*yNuGrid.max())
plt.show()



##================================    convolve    ==================================
convY = np.convolve(yNuGrid,gaussY,'same')
plt.plot(xNuGrid,convY)
plt.show()



##==================   Convert back to wavelength space   ==========================
xWvl = c/xNuGrid*1E8
yWvl = convY/(xWvl**2)*3E-5 #convert Fnu(Jy) to Flambda
#yWvl *= fluxConservationFactor #renormalize so flux is conserved
plt.plot(xWvl[xWvl<16000], yWvl[xWvl<16000])
plt.plot(x,y/max(y)*max(yWvl))
plt.show()
'''

