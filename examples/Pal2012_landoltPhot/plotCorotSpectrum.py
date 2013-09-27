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
        filty = filty/max(filty)#*0.8218 #normalize to our Johnson R filter max Transmission
    if fname == 'Vfilter.txt':
        filty = filty/max(filty)#*0.8623 #normalize to our Johnson V filter max transmission

    if (filtx[-1] < filtx[0]): #if array goes high to low, reverse the order to the integration does not get a negative
        filtx = filtx[::-1]
        filty = filty[::-1]

    filtWidth = filtx[-1]-filtx[0]
    print "filter width = ",filtWidth
    filtCorrection = integrate.simps(filty,x=filtx)/filtWidth #calculate correction factor for filter width
    print "filter correction = ", filtCorrection
    return filtx, filty, filtWidth, filtCorrection

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

#LOAD A SINGLE FILE
#FileName = '/home/srmeeker/scratch/standards/crab_fit_11.npz'
#FileName = '/home/srmeeker/scratch/standards/sdssj0926/sdssj0926_fit_6.npz'
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

diam = 500 #5 meter telescope
area = np.pi * ((diam/2.0)**2 -(100)**2) #Palomar secondary is ~1m radius from Serabyn 2007
curve/= area #spectrum is now in counts/s/Angs/cm^2



#SETUP PLOTTING
fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlim(3500,13000)
#plt.ylim(0,0.001)

fname = 'throughput.txt'
fdata = np.loadtxt(fname,dtype=float)
x = np.array(fdata[:,0])
tp = np.array(fdata[:,1])

ax.plot(wvls,curve)
#ax.plot(x,y)
ax.set_title('Absolute Spectrum of  '+FileName.split('/')[-1].split('_')[0])
#ax.set_yscale('log')

#add points for crab photometry
# bands are [U,B,V,R,I,J]
#mags = [16.68,17.22,16.64,16.14,15.61,14.72] #from Sandberg and Sollerman, 2009, Crab+Knot mags

#errors = [0.03,0.02,0.02,0.01,0.01,0.03]

centers = [3600,4400,5500,6400,7900,12500]
zeros = [759,1461,999,726,487,194]

#corot Mags
mags = [15,15.79,15.00,14.472,14.051,13.441]
errors = [0,0,0.1,0.048,0.03,0.024]

#plt.show()

#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(wvls,curve/tp)
#ax.set_ylim(min((curve/tp)[wvls>3800])/2.0,max((curve/tp)[wvls>3800])*1.5)
#ax.set_xlim(3800,12500)

#Plot COROT18 spectrum for comparison
#load v band filter transmission
fname = 'Vfilter.txt'
tx, ty, filtWidth, filtCorrection = loadFilter(fname)

'''
#load synth spectrum
fname = 'G9V.spec'
fdata = np.loadtxt(fname,dtype=float)
x = np.array(fdata[:,0])
y = np.array(fdata[:,1]) #flux in ergs/s/cm^2/Angs
y = y/(h*c/x) #flux in counts/s/cm^2/Angs
plt.plot(x,y,color='red',alpha=0.4)
synthx = x
synthSpectrum = y
#get v-band portion of synthetic spectrum
transmission = interpolate.griddata(tx,ty,x)
print transmission
transmission[np.isnan(transmission)] = 0.0
print transmission
vcounts = y*transmission
plt.plot(x,vcounts,color='black')
print vcounts
totalSynthVCounts = integrate.simps(vcounts,x=x) #gives total counts in counts/s/cm^2 after integrating over Angs
filtWidth = tx[-1]-tx[0]
print "filter width = ",filtWidth
filtCorrection = integrate.simps(transmission,x=x)/filtWidth #calculate correction factor for filter width
print "filter correction = ", filtCorrection
totalSynthVCounts = totalSynthVCounts/(filtWidth)/filtCorrection #divide by width of filter to get per Angstrom
#corot parameters
'''
corotVMag = 14.99
corotVCounts = (10.0**(corotVMag/(-2.5)))*zeros[2]#*filtCorrection #still in counts/s/cm^2/Angs

corotMagCounts = (10.0**(np.array(mags)/(-2.5)))*np.array(zeros)
counterrors = corotMagCounts-10**((np.array(mags)+np.array(errors))/(-2.5))*np.array(zeros)

'''
corotScale = corotVCounts/totalSynthVCounts
print totalSynthVCounts
print corotVCounts
print "corot Scaling = ",corotScale
'''


fig = plt.figure()
ax = fig.add_subplot(111)
#plt.plot(wvls,curve/tp)
#plt.errorbar(wvls,curve/tp,np.sqrt(curve*300*area*binWidths)/300/area/binWidths/tp)
#plt.step(wvls,curve/tp,color='black',where='mid')
plt.step(wvlBinEdges[1:],curve/tp, color='black', where='pre')
#plt.plot(synthx,synthSpectrum*corotScale,color="red",alpha = 0.4)
#bin corot spectrum to our energy resolution
energyBinWidth = 0.3936
wvlStart = 3000
wvlStop = 13000
wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth,wvlStart,wvlStop)
#rebinned = rebin(synthx,(synthSpectrum*corotScale),wvlBinEdges)
#plt.plot(rebinned[:,0],rebinned[:,1],color="black")

plt.errorbar(centers,corotMagCounts,yerr=counterrors,fmt='o')

ax.set_ylim(min((curve/tp)[wvls>3800])/2.0,max((curve/tp)[wvls>3800])*1.5)
ax.set_xlim(3800,13000)

ax.set_title("Corot 18")
ax.set_xlabel("Wavelength [$\AA$]")
ax.set_ylabel("Flux [photons/s/cm$^{2}$/$\AA$]")
plt.legend(["ARCONS Spectrum","BVRIJ Magnitudes"],numpoints=1)

#ax.set_yscale("log")
plt.show()

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

