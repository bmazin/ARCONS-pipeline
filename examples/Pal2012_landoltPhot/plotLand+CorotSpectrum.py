from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from util.rebin import rebin
from matplotlib import rcParams
#from matplotlib.backends.backend_pdf import PdfPages
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


if len(sys.argv) >1:
    fileNum = str(sys.argv[1])
else:
    fileNum = '0'

# common setup for matplotlib
params = {'axes.labelsize': 14,
          'text.fontsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 12,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'xtick.minor.pad': 6,
          'ytick.labelsize': 12,
          'lines.linewidth':1,
          'lines.markersize'  : 6,
          # free font similar to Helvetica
          'font.family':'FreeSans'}

rcParams.update(params)

c=3.00E18 #Angs/s
h=6.626E-27 #erg*s

#LOAD ALL SDSSj0926/CRAB FILES

NumFrames = 31
#nFiles = 13 #j0926
#nFiles=25 #crab
nFiles = 1 #try with only 1 file
curves = np.zeros((nFiles,NumFrames),dtype=float)

for k in xrange(nFiles):
    FileName = '/home/srmeeker/scratch/standards/Landolt9542_fit_%s.npz'%(fileNum)
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

    diam = 510.55 #5 meter telescope
    area = np.pi * ((diam/2.0)**2 -(183/2.0)**2) #secondary obstruction diameter 1.83m
    print k
    print curve
    curves[k] = curve
print curves
finalCurve = np.zeros((NumFrames),dtype=float)
for j in xrange(NumFrames):
    finalCurve[j] = np.median(curves[:,j])

curve = finalCurve


#LOAD A SINGLE FILE
#FileName = '/home/srmeeker/scratch/standards/crabNight2_fit_1.npz'
#FileName = '/home/srmeeker/scratch/standards/sdssj0926/sdssj0926_fit_6.npz'
#FileName = '/home/srmeeker/scratch/standards/corot18/corot18_fit.npz'
NumFrames = 31

t = np.load(FileName)
#energyBinWidth = 0.3936
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

curve = 2*np.pi*amps*(widths**2) #spectrum of oberved object in counts/s
curve /= binWidths #spectrum is now in counts/s/Angs

diam = 510.55 #5 meter telescope
area = np.pi * ((diam/2.0)**2 -(183/2.0)**2) #secondary obstruction diameter 1.83m
curve/= area #spectrum is now in counts/s/Angs/cm^2


#SETUP PLOTTING
#fig = plt.figure()
#ax = fig.add_subplot(111)
#plt.xlim(3500,13000)
#plt.ylim(0,0.001)

#SETUP EPS PLOTTING
#matplotlib.rcParams.update({'font.size':12, 'font.family': 'sans-serif','sans-serif':['Helvetica']})
plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('text',usetex=True)

fname = 'throughput_high.txt'
fdata = np.loadtxt(fname,dtype=float)
xh = np.array(fdata[:,0])
tph = np.array(fdata[:,1])

fname = 'throughput.txt'
fdata = np.loadtxt(fname,dtype=float)
x = np.array(fdata[:,0])
tp = np.array(fdata[:,1])


#add points for crab photometry
# bands are [U,B,V,R,I,J]
mags = [14.28,15.391,15.606,15.725,15.905,16.130] #noao.edu
centers = [3600,4400,5500,6400,7900,12500]
errors = [0.0046,0.0044,0.0058,0.0048,0.028,0.066]
zeros = [759,1461,999,726,487,194]

#centers = [3600,4400,5500,7000,9000,12500]
#zeros = [4.35e-09,7.20e-09,3.92e-09,1.76e-09,8.30e-10,3.40e-10] #ergs/s/cm^2/Angs

vx,vy,vwidth,vcorr = loadFilter('Vfilter.txt')
rx,ry,rwidth,rcorr = loadFilter('Rfilter.txt')
corrections = np.array([1,1,vcorr,rcorr,1,1])
#corrections = np.array([1,1,1,1,1,1])

counts = 10**(np.array(mags)/(-2.5))*np.array(zeros) #* corrections

counterrors = counts-10**((np.array(mags)-np.array(errors))/(-2.5))*np.array(zeros)#*corrections
print "Calculated Landolt counts from photometry = "
print counts
print counterrors
percents = counterrors/counts * 100
print "percents (errors/counts * 100) = ", percents

#plt.show()
plotDir = "/home/srmeeker/scratch/standards"
#plotFileName = "Land+CorotSpec.pdf"
plotFileName = "Land+CorotSpec.eps"
fullFluxPlotFileName = os.path.join(plotDir,plotFileName)

#pp = PdfPages(fullFluxPlotFileName)

fig = plt.figure()
ax = plt.subplot(211)
plt.subplots_adjust(hspace=0.01)

#put into ergs
curve = curve*(h*c/wvls)*1e15

plt.step(wvlBinEdges[1:],curve/tp, color='black', where='pre')

ax.set_ylim(min((curve/tp)[wvls>3800])/5.0,max((curve/tp)[wvls>3800])*1.0)
#ax.set_ylim(min((curve/tp)[wvls>3800])/3,max((curve/tp)[wvls>3800])/2.0)
ax.set_xlim(4000,11000)

#ax.set_xscale('log')

ax.yaxis.labelpad = 15
ax.xaxis.labelpad = -1
ax.xaxis.set_visible(False)

#put crab photometry points on plot
counts *= (h*c/np.array(centers))*1e15
counterrors *= (h*c/np.array(centers))*1e15

plt.errorbar(centers,counts,yerr=counterrors,fmt='o')

#ax.set_title("Crab Pulsar")
ax.set_xlabel("Wavelength [$\AA$]")
#ax.set_ylabel("Flux [photons/s/cm$^{2}$/$\AA$]")
plt.legend(["ARCONS Spectrum of Landolt 95-42","BVRI Fluxes"],numpoints=1)


FileName = '/home/srmeeker/scratch/standards/corot18_fit_0.npz'
NumFrames = 31

t = np.load(FileName)
#energyBinWidth = 0.3936
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

curve = 2*np.pi*amps*(widths**2) #spectrum of oberved object in counts/s
curve /= binWidths #spectrum is now in counts/s/Angs
curve/= area #spectrum is now in counts/s/Angs/cm^2

# bands are [U,B,V,R,I,J]
#corot Mags
mags = [15,15.79,15.00,14.472,14.051,13.441]
errors = [0,0,0.1,0.048,0.03,0.024]

#Plot COROT18 spectrum for comparison
#load v band filter transmission
fname = 'Vfilter.txt'
tx, ty, filtWidth, filtCorrection = loadFilter(fname)

corotVMag = 14.99
corotVCounts = (10.0**(corotVMag/(-2.5)))*zeros[2]#*filtCorrection #still in counts/s/cm^2/Angs

corotMagCounts = (10.0**(np.array(mags)/(-2.5)))*np.array(zeros)
counterrors = corotMagCounts-10**((np.array(mags)+np.array(errors))/(-2.5))*np.array(zeros)

ax = fig.add_subplot(212)

#put into ergs
curve = curve*(h*c/wvls)*1e15

plt.step(wvlBinEdges[1:],curve/tp, color='black', where='pre')

corotMagCounts *= (h*c/np.array(centers))*1e15
counterrors *= (h*c/np.array(centers))*1e15

plt.errorbar(centers,corotMagCounts,yerr=counterrors,fmt='o')

ax.set_ylim(min((curve/tp)[wvls>3800])/2.0,max((curve/tp)[wvls>3800])*1.4)
ax.set_xlim(4000,11000)
#ax.set_title("Corot 18")
plt.legend(["ARCONS Spectrum of Corot 18","BVRI Fluxes"],numpoints=1)

#ax.set_xscale("log")

xl = plt.xlabel(ur"Wavelength [$\AA$]")
yl = plt.ylabel("F$_\lambda$ [$10^{-15}$ ergs/sec/cm$^{2}$/$\AA$]")
yl.set_y(1)
ax.yaxis.labelpad = 15
ax.xaxis.labelpad = 10
ax.tick_params(axis='x',which='minor',bottom='on',labelbottom='on')
ax.xaxis.set_minor_formatter(plt.FormatStrFormatter('%d'))
ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))

#plt.show()

plt.savefig(fullFluxPlotFileName,format='eps')

#pp.savefig()
#pp.close()

'''
#output spectrum
outspec = synthSpectrum*corotScale
outx = np.array(range(20000))
print outx
outspec = interpolate.griddata(synthx,outspec,outx)
outspec[np.isnan(outspec)]=0
outspec = outspec*(h*c/outx)
print outspec
fig = plt.figure()
ax = fig.add_subplot(111)
plt.plot(outx,outspec)
ax.set_xlim(0,20000)
plt.show()

outarr = np.empty((len(outx),1),dtype=float)
outarr[:,0]=outspec
np.savetxt("G9V_spec_rebinned.txt", outarr)
'''






















