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
#nFiles=25 #only files in bad position
#nFiles = 34 #try with only 1 file
#offset=0
nFiles = 9 #only 9 from good position
offset = 25
curves = np.zeros((nFiles,NumFrames),dtype=float)

for k in xrange(nFiles):
    FileName = '/home/srmeeker/scratch/standards/crabNight2_fit_%s.npz'%(k+offset)
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

diam = 500 #5 meter telescope
area = np.pi * ((diam/2.0)**2 -(100)**2) #Palomar secondary is ~1m radius from Serabyn 2007
curve/= area #spectrum is now in counts/s/Angs/cm^2


#SETUP PLOTTING
#fig = plt.figure()
#ax = fig.add_subplot(111)
#plt.xlim(3500,13000)
#plt.ylim(0,0.001)

fname = 'throughput_high.txt'
fdata = np.loadtxt(fname,dtype=float)
xh = np.array(fdata[:,0])
tph = np.array(fdata[:,1])

fname = 'throughput.txt'
fdata = np.loadtxt(fname,dtype=float)
x = np.array(fdata[:,0])
tp = np.array(fdata[:,1])

#ax.plot(wvls,curve)
#ax.plot(x,y)
#ax.set_title('Absolute Spectrum of  '+FileName.split('/')[-1].split('_')[0])
#ax.set_yscale('log')

#add points for crab photometry
# bands are [U,B,V,R,I,J]
mags = [16.68,17.22,16.64,16.14,15.61,14.72] #from Sandberg and Sollerman, 2009, Crab+Knot mags
centers = [3600,4400,5500,6400,7900,12500]
errors = [0.03,0.02,0.02,0.01,0.01,0.03]
zeros = [759,1461,999,726,487,194]

vx,vy,vwidth,vcorr = loadFilter('Vfilter.txt')
rx,ry,rwidth,rcorr = loadFilter('Rfilter.txt')
corrections = np.array([1,1,vcorr,rcorr,1,1])
#corrections = np.array([1,1,1,1,1,1])

counts = 10**(np.array(mags)/(-2.5))*np.array(zeros) #* corrections

#USING BEN'S MAGNITUDE CALC CODE
#janZeros = [1810,4260,3640,3080,2550,1600]
#dLambdaOverLambda = [0.15,0.22,0.16,0.23,0.19,0.16]
#lambda_c = [0.36,0.44,0.55,0.64,0.79,1.26]
#counts = np.array(janZeros) * 10.0**(-1*np.array(mags)/2.5)*1.517e7*np.array(dLambdaOverLambda)/10000.0 #photons/sec/cm^2

counterrors = counts-10**((np.array(mags)-np.array(errors))/(-2.5))*np.array(zeros)#*corrections
print "Calculated Crab counts from photometry = "
print counts
print counterrors
percents = counterrors/counts * 100
print "percents (errors/counts * 100) = ", percents

#plt.show()
#plotDir = "/home/srmeeker/scratch/standards"
#plotFileName = "crabNight2_spec_%s.pdf"%(fileNum)
#fullFluxPlotFileName = os.path.join(plotDir,plotFileName)
#pp = PdfPages(fullFluxPlotFileName)
#matplotlib.rcParams['font.size']=6

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(wvls,curve/tp)
#ax.plot(wvls,curve/tph,'grey')
plt.step(wvlBinEdges[1:],curve/tp, color='black', where='pre')
ax.set_ylim(min((curve/tp)[wvls>3800])/2.0,max((curve/tp)[wvls>3800])*1.5)
#ax.set_ylim(min((curve/tp)[wvls>3800])/3,max((curve/tp)[wvls>3800])/2.0)
ax.set_xlim(3800,12500)

#put crab photometry points on plot
plt.errorbar(centers,counts,yerr=counterrors,fmt='o')

ax.set_title("Crab Pulsar")
ax.set_xlabel("Wavelength [$\AA$]")
ax.set_ylabel("Flux [photons/s/cm$^{2}$/$\AA$]")
plt.legend(["ARCONS Spectrum","BVRIJ Magnitudes"],numpoints=1)

#ax.set_yscale("log")
plt.show()

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






















