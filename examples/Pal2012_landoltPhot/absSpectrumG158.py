from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from headers import pipelineFlags
from util.rebin import rebin
import matplotlib
from mpltools	import style
from scipy import interpolate
from scipy.optimize.minpack import curve_fit
from math import exp
import figureHeader

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


c=3.00E18 #Angs/s
h=6.626E-27 #erg*s

FileName = '/home/srmeeker/scratch/standards/G158-100_fit_0.npz'
NumFrames = 31

t = np.load(FileName)

energyBinWidth = 0.1
wvlStart = 3000
wvlStop = 13000
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
curve /= binWidths #spectrum is now in counts/s/Angs

diam = 510.55 #5 meter telescope
area = np.pi * ((diam/2.0)**2 -(183/2.0)**2) #secondary obstruction diameter 1.83m
curve/= area #spectrum is now in counts/s/Angs/cm^2

#SETUP PLOTTING
fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlim(3700,16000)
#plt.ylim(0,0.001)

#LOAD KNOWN SPECTRUM OF STANDARD
objectName = "G158-100"
#import the known spectrum of the calibrator and rebin to the histogram parameters given
#must be imported into array with dtype float so division later does not have error
#std = MKIDStd.MKIDStd()
#a = std.load(objectName)
#x = a[:,0]
#y = np.array(a[:,1]) #std object spectrum in ergs/s/Angs/cm^2
#print y
#y*=(1E-16) #G158 saved in units of 10^-16 ergs/s/Angs/cm^2
#print y
#convert from ergs/s to counts/s

#Manually load G158 standard file. MKIDStd units may be wrong.
#scale = (1.0E16) #for .dat files strange unit scaling
#fname = 'fg158_100.dat'# flux in 10^16 ergs/s/cm^2/Angs

fname = "G158_Filippenko_1984.txt" #flux in micro-janksys
scale = 10.0**6.0 #convert to janskys

fdata = np.loadtxt(fname,dtype=float)
x = np.array(fdata[:,0]) #angstroms
#y = np.array(fdata[:,1])/scale #flux in ergs/s/cm^2/Angs

y = np.array(fdata[:,1])/scale #convert to flux in Janskys
y = y*(3e-5)/(x**2) #flux in ergs/s/cm^2/Angs

print y
y = y/(h*c/x) #flux in counts/s/cm^2/Angs
print y
print "Loaded standard spectrum of %s"%(objectName)

newwl, newflux = cleanSpectrum(x,y,objectName,wvlBinEdges)
print newwl
print newflux
plt.plot(newwl,newflux)
plt.plot(x,y)
plt.title("Standard spectrum with exponential tail fit")
plt.legend(['Exponential Tail','Rebinned Spectrum','Final Refitted spectrum','Original Spectrum'],'br')
plt.show()

newa = rebin(newwl,newflux,wvlBinEdges)
x = newa[:,0]
y = newa[:,1]

#print x
#print y
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wvls,curve)
ax.plot(x,y)
ax.set_title('Absolute Spectrum of  '+FileName.split('/')[-1].split('_')[0])
ax.set_yscale('log')
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(wvls,curve/y)
ax.set_ylim(0,0.02)
ax.set_xlim(3750,12500)
ax.set_title('Throughput')
plt.show()

# try new plotting style for pipeline paper
#style.use('ggplot')
fig = plt.figure()
ax = fig.add_subplot(111)
#Plot 2012 measured throughput data inverse: flux cal correction curve
ax.plot(wvls,1/(curve/y),'black',linewidth=3)#, label='On-sky Total Throughput (2012)')
ax.set_xlim(4000,11000)
ax.set_ylim(0,250)
plt.xlabel(ur"Wavelength [$\AA$]")
plt.ylabel(ur"Spectral Calibration Curve")
plt.savefig("FluxCal_SensitivityCurve.eps",format='eps')
outdir = '/home/srmeeker/scratch/junk/'
np.savez(outdir+'%s_throughput.npz'%(objectName.strip()),throughput=curve/y,wvls=wvls)

#MAKE FLUXSOL FILE WITH ABSOLUTE FLUX CALIBRATION
'''
fluxCalFileName = "fluxsol_absolute_021727.h5"
if os.path.isabs(fluxCalFileName) == True:
    fullFluxCalFileName = fluxCalFileName
else:
    scratchDir = os.getenv('INTERM_PATH')
    fluxDir = os.path.join(scratchDir,'fluxCalSolnFiles/PAL2012/20121211/')
    fullFluxCalFileName = os.path.join(fluxDir,fluxCalFileName)

try:
    fluxCalFile = tables.openFile(fullFluxCalFileName,mode='w')
except:
    print 'Error: Couldn\'t create flux cal file, ',fullFluxCalFileName

fluxFactors = y/curve
fluxFlags = np.empty(np.shape(fluxFactors),dtype='int')
fluxFlags.fill(pipelineFlags.fluxCal['good'])   #Initialise flag array filled with 'good' flags. JvE 5/1/2013.
#set factors that will cause trouble to 1
#self.fluxFlags[self.fluxFactors == np.inf] = 1
fluxFlags[fluxFactors == np.inf] = pipelineFlags.fluxCal['infWeight']   #Modified to use flag dictionary - JvE 5/1/2013
fluxFactors[fluxFactors == np.inf]=1.0
fluxFactors[np.isnan(fluxFactors)]=1.0
fluxFlags[np.isnan(fluxFactors)] = pipelineFlags.fluxCal['nanWeight']   #Modified to use flag dictionary - JvE 5/1/2013
fluxFlags[fluxFactors <= 0]=pipelineFlags.fluxCal['LEzeroWeight']   #Modified to use flag dictionary - JvE 5/1/2013
fluxFactors[fluxFactors <= 0]=1.0

calgroup = fluxCalFile.createGroup(fluxCalFile.root,'fluxcal','Table of flux calibration weights by wavelength')
caltable = tables.Array(calgroup,'weights',object=fluxFactors,title='Flux calibration Weights indexed by wavelengthBin')
flagtable = tables.Array(calgroup,'flags',object=fluxFlags,title='Flux cal flags indexed by wavelengthBin. 0 is Good')
bintable = tables.Array(calgroup,'wavelengthBins',object=wvlBinEdges,title='Wavelength bin edges corresponding to third dimension of weights array')
fluxCalFile.flush()
fluxCalFile.close()
print "Finished Flux Cal, written to %s"%(fullFluxCalFileName)
'''

