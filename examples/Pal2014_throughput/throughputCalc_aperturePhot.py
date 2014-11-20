from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util.ObsFile import ObsFile
from util import MKIDStd
from util.readDict import readDict
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
        #weight = [0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.7]
        weight = 0.7*np.ones((len(re_wl)))
        #print len(weight)
        #weight = re_flux/min(re_flux)
        #weight = 1.0/weight
        #weight = weight/max(weight)
        print weight
        f = interpolate.splrep(re_wl,re_flux,w=weight,k=2,s=0)#max(re_flux)**300)
        new_flux = interpolate.splev(new_wl,f,der=0)
        return new_wl, new_flux

def aperture(startpx,startpy,radius=7):
    r = radius
    length = 2*r 
    height = length
    allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
    ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
    pixx = []
    pixy = []
    mask=np.ones((46,44))
    for x in allx:
        for y in ally:
            if (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 <= (r)**2 and 0 <= y and y < 46 and 0 <= x and x < 44:
                mask[y,x]=0.
    return mask


c=3.00E10 #cm/s
h=6.626E-27 #erg*s
k=1.3806488E-16 #erg/K

param = readDict()
if len(sys.argv)<2:
    print "Provide file name to fit.  Syntax >>python throughputCalc_aperturePhot.py objectparams.dict [filenumber]"
    sys.exit(1)

#read in parameter file as command line argument
param.read_from_file(sys.argv[1])

#provide optional file number if the object in the param file has alternate .npz files to be specified individually
fileNum = None
if len(sys.argv)>2:
    fileNum = "_"+str(sys.argv[2])

npzfitpsf = param['npzfitpsf']

if fileNum != None:
    npzfitpsf = npzfitpsf.split('.')[0]+fileNum+'.'+npzfitpsf.split('.')[1]

FramesPerFile = param['FramesPerFile']
guessX = param['guessX'][0]
guessY = param['guessY'][0]
objectName = param['objectName']

#SETUP PLOTTING
#matplotlib.rcParams.update({'font.size':12, 'font.family': 'sans-serif','sans-serif':['Helvetica']})
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('text',usetex=True)

t = np.load(npzfitpsf)

energyBinWidth = 0.1
wvlStart = 3000
wvlStop = 25000
wvlBinEdges = ObsFile.makeWvlBins(energyBinWidth,wvlStart,wvlStop)
nWvlBins = len(wvlBinEdges)-1
binWidths = np.empty(nWvlBins)
for i in xrange(nWvlBins):
    binWidths[i] = wvlBinEdges[i+1]-wvlBinEdges[i]
#print binWidths
cube = t['cube']
wvls = t['wvls']
curve=[]

for iFrame in range(0,np.shape(cube)[0]):
    print wvls[iFrame]
    frame = cube[iFrame,:,:]
    objectLight = np.zeros((np.shape(frame)),dtype = float)
    skyLight = np.zeros((np.shape(frame)),dtype = float)

    nanMask = np.isnan(frame)
    err = np.sqrt(frame)

    apertureMask = aperture(guessX,guessY,radius=14)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(wvls[iFrame])
    im = ax.matshow(apertureMask, cmap = cm.get_cmap('rainbow'))
    fig.colorbar(im)
    #plt.show()

    objectLight[:,:] = frame[:,:]
    objectLight[apertureMask!=0] = np.nan
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(wvls[iFrame])
    im = ax.matshow(objectLight, cmap = cm.get_cmap('rainbow'))
    fig.colorbar(im)
    #plt.show()

    skyLight[:,:] = frame[:,:]
    skyLight[apertureMask==0] = np.nan
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(wvls[iFrame])
    im = ax.matshow(skyLight, cmap = cm.get_cmap('rainbow'))
    fig.colorbar(im)
    #plt.show()

    skyLightm = np.ma.masked_array(skyLight, [np.isnan(x) for x in skyLight])
    skyLevel = np.ma.median(skyLightm)
    print "sky median = ", skyLevel
    skySubtracted = objectLight - skyLevel
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(wvls[iFrame])
    im = ax.matshow(skySubtracted, cmap = cm.get_cmap('rainbow'))
    fig.colorbar(im)
    #plt.show()

    countsPerSec = np.nansum(skySubtracted)
    print "aperture sum = ", countsPerSec
    curve.append(countsPerSec)

curve = np.array(curve)
curve /= binWidths #spectrum is now in counts/s/Angs

np.savez('%s_CountsPerSecondPerBin.npz'%(objectName.strip()),counts=curve,binEdges=wvlBinEdges,binWidths=binWidths,centerWvls=wvls)

diam = 510.55 #5 meter telescope
area = np.pi * ((diam/2.0)**2 -(183/2.0)**2) #secondary obstruction diameter 1.83m
curve/= area #spectrum is now in counts/s/Angs/cm^2

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(4000,25000)
#ax.set_ylim(0,1E-5)
plt.plot(wvls, curve)
plt.xlabel(ur"Wavelength [\AA]")
plt.ylabel(ur"ARCONS measured Spectrum (Counts/s/cm2/Angs)")
plt.savefig("FluxCal_RawSpectrumCounts.eps",format='eps')

#convert ARCONS measured spectrum to ergs/s/cm^2/A
curve*=h*(c*1E8)/wvls

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(4000,25000)
#ax.set_ylim(0,0.35E-16)
plt.plot(wvls, curve)
plt.xlabel(ur"Wavelength [\AA]")
plt.ylabel(ur"ARCONS measured Spectrum (Flambda)")
plt.savefig("FluxCal_RawSpectrumFlux.eps",format='eps')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(4000,25000)
#Begin MKIDStd loading
#objectName = "lds749b"
#import the known spectrum of the calibrator and rebin to the histogram parameters given
#must be imported into array with dtype float so division later does not have error
std = MKIDStd.MKIDStd()
a = std.load(objectName)
a = std.countsToErgs(a)
x = a[:,0]
y = np.array(a[:,1]) #std object spectrum in counts/s/Angs/cm^2
#convert from counts/s to ergs/s
#y=y*(h*(c*1E8)/x)
#End MKIDStd loading

plt.plot(x,y*1E15,linewidth=1,color='grey',alpha=0.75)
newwl, newflux = cleanSpectrum(x,y,objectName,wvlBinEdges)
print newwl
print newflux
#plt.plot(newwl,newflux*1E15,color = 'red')
#plt.show()

newa = rebin(newwl,newflux,wvlBinEdges)

x = newa[:,0]
y = newa[:,1]

plt.step(x,y*1E15,color = 'black',where='mid')
plt.legend(['LDS749B Spectrum','BB Fit','Rebinned Std Spectrum','Resampled Std Spectrum'],'upper right', numpoints=1)
plt.xlabel(ur"Wavelength [\AA]")
plt.ylabel(ur"Flux [10$^{-15}$ ergs/s/cm$^{2}$/\AA]")
plt.ylim(0,15)
plt.savefig('FluxCal_StdSpectrum.eps',format='eps')
#plt.show()

bvrwvls = [4450, 5510, 6580]#center wvls for b v and r Johnson filters
widths = [940/2.0, 880/2.0,1380/2.0] #filter FWHMs
bvrthru = [.23,.24,.35] #as calculated by Pal2013 throughput code
errors = [.04,.04,.04] #as calculated by Pal2013 throughput code (need to get real numbers, these are fake now)

#load QE file for plotting
QEFileName = "avgQE_20131125-202654.txt"
QEfile = os.environ['ARCONS_PIPELINE_PATH']+'/util/data/'+QEFileName
fdata = np.loadtxt(QEfile,dtype=float)
qewvls = np.array(fdata[:,0])*10.0 #convert from nm to Angstroms
QEcurve = np.array(fdata[:,1])


fig = plt.figure()
ax = fig.add_subplot(111)
#plt.errorbar(bvrwvls, np.array(bvrthru)*100, xerr = widths, yerr=np.array(errors)*100, fmt='o',color='black')
#ax.plot(qewvls, QEcurve*100,linestyle="--",color='black')
#ax.plot(wvls,(curve/y)*100,'black')
#ax.set_ylim(4E-3,0.04)
ax.set_ylim(0,5)
ax.set_xlim(4000,11000)
plt.plot(wvls, (curve/y)*100)
plt.xlabel(ur"Wavelength [\AA]")
plt.ylabel(ur"Throughput [\%]")

multqewvls = [4500, 5500, 6500]
multqe = [QEcurve[qewvls==4500][0]*bvrthru[0], QEcurve[qewvls==5500][0]*bvrthru[1], QEcurve[qewvls==6500][0]*bvrthru[2]]
#plt.errorbar(multqewvls, np.array(multqe)*100, color='blue',fmt='o')
print multqe

#plt.title("Sensitivity Curve")
#plt.legend(['Telescope BVR Throughput','ARCONS QE', 'Total QE (Telescope + ARCONS)'],'upper right', numpoints=1)
#ax.set_yscale('log')
plt.savefig("FluxCal_SensitivityCurve.eps",format='eps')
#plt.show()

np.savez('%s_throughput.npz'%(objectName.strip()),throughput=curve/y,wvls=wvls)

