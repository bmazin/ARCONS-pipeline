# -*- coding: utf-8 -*-
from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from util import MKIDStd
from util.rebin import rebin
import matplotlib
from mpltools	import style
from scipy import interpolate
from scipy.optimize.minpack import curve_fit
from matplotlib.ticker import ScalarFormatter 
from numpy import exp
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
        #R=E/dE -> dE = E/R
        dE = 0.3936 #eV
        start = 1000 #Angs
        stop = 25000 #Angs
        enBins = ObsFile.makeWvlBins(dE,start,stop)
        rebinned = rebin(wl,flux,enBins)
        re_wl = rebinned[:,0]
        re_flux = rebinned[:,1]
        plt.plot(re_wl,re_flux*1E15,linestyle="o", marker="o",markersize=10, color="blue") #plot rebinned spectrum with exp tail
        
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

# try new plotting style for pipeline paper
#style.use('ggplot')

c=3.00E10 #cm/s
h=6.626E-27 #erg*s
k=1.3806488E-16 #erg/K

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
#matplotlib.rcParams.update({'font.size':12, 'font.family': 'sans-serif','sans-serif':['Helvetica']})

#ignore these parameters for Julian's figureHeader params
#plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
#plt.rc('text',usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111)
plt.xlim(4000,11000)
#plt.ylim(0,0.001)

#Begin MKIDStd loading
objectName = "G158-100"
#import the known spectrum of the calibrator and rebin to the histogram parameters given
#must be imported into array with dtype float so division later does not have error
std = MKIDStd.MKIDStd()
a = std.load(objectName)
x = a[:,0]
y = np.array(a[:,1]) #std object spectrum in counts/s/Angs/cm^2
#convert from counts/s to ergs/s
y=y*(h*(c*1E8)/x)
#End MKIDStd loading

'''#Begin manual standard loading block
#Manually load G158 standard file. MKIDStd units may be wrong.
#scale = (1.0E16) #for .dat files strange unit scaling
#fname = 'fg158_100.dat'# flux in 10^16 ergs/s/cm^2/Angs

fname = "G158_Filippenko_1984.txt" #flux in micro-janksys
scale = 1.0E6 #for converting to janskys from microjanskys

fdata = np.loadtxt(fname,dtype=float)
x = np.array(fdata[:,0]) #angstroms
#y = np.array(fdata[:,1])/scale #flux in ergs/s/cm^2/Angs

y = np.array(fdata[:,1])/scale #convert to flux in Janskys
y = y*(3E-5)/(x**2) #flux in ergs/s/cm^2/Angs

#y = y/(h*c/x) #flux in counts/s/cm^2/Angs
#print y
print "Loaded standard spectrum of %s"%(objectName)
'''#End manual standard loading block

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

#convert ARCONS measured spectrum to ergs/s/cm^2/A
curve*=h*(c*1E8)/x

#print x
#print y

#fig = plt.figure()
#ax = fig.add_subplot(111)
#plt.plot(wvls,curve*1E15) #plot ARCONS measured spectrum
#ax.plot(x,y)
#ax.set_title('Absolute Spectrum of  '+FileName.split('/')[-1].split('_')[0])
#ax.set_yscale('log')

#plt.title("Spectrophotometric Standard Spectrum")
#plt.legend(['Original Std Spectrum','BB Fit','Rebinned Std Spectrum','Resampled Std Spectrum','ARCONS Measured Spectrum'],'right', numpoints=1)
plt.legend(['G158-100 Spectrum','BB Fit','Rebinned Std Spectrum','Resampled Std Spectrum'],'upper right', numpoints=1)
plt.xlabel(ur"Wavelength (\r{A})")
plt.ylabel(ur"Flux (10$^{-15}$ ergs s$^{-1}$ cm$^{-2}$ \r{A}$^{-1}$)")
plt.ylim(1,5)
plt.savefig('FluxCal_StdSpectrum.eps',format='eps')
#plt.show()

bvrwvls = [4450, 5510, 6580]#center wvls for b v and r Johnson filters
widths = [940/2.0, 880/2.0,1380/2.0] #filter FWHMs
#bvrthru = [.237,.258,.35] #as calculated by Pal2013 throughput code
#Readjusted standard star intensity by factor of 1.2 to get it to match known magnitudes better when multiplied through our filters
bvrthru = [.199, .215, .273]
errors = [.03,.05,.045] #as calculated by Pal2013 throughput code

#load QE file for plotting
QEFileName = "avgQE_20131125-202654.txt"
QEfile = os.environ['ARCONS_PIPELINE_PATH']+'/util/data/'+QEFileName
fdata = np.loadtxt(QEfile,dtype=float)
qewvls = np.array(fdata[:,0])*10.0 #convert from nm to Angstroms
QEcurve = np.array(fdata[:,1])

fig = plt.figure()
ax = fig.add_subplot(111)

####Plot telescope throughput as measured with photo diode
plt.errorbar(bvrwvls, np.array(bvrthru)*100, xerr = widths, yerr=np.array(errors)*100, fmt='o',color='black',markersize=8)#,label='Telescope BVR Throughput (2013)')

####Plot ARCONS LAB QE 2013
ax.plot(qewvls, QEcurve*100,linestyle="--",color='black')#, label='ARCONS QE (2013)')
plt.fill_between(qewvls,QEcurve*100+1,QEcurve*100-1,color='LightSteelBlue',alpha=0.2)

####plot calculated total throughput as telescope tp * arcons lab QE
multqewvls = [4500, 5500, 6500]
QEcurvePoints = [QEcurve[qewvls==4500][0], QEcurve[qewvls==5500][0], QEcurve[qewvls==6500][0]]
multqe = np.array(QEcurvePoints) * np.array(bvrthru)

#estimate errors for QE data just as 1% for now. Don't save that info in the measurement, but QE alignment could make up 2% alone.
#sig^2 = df/da^2*sig_a^2 where f is QEcurvePoints*bvrthru and sig_a is errors
multqeErrs = np.sqrt(np.array(errors)**2*np.array(QEcurvePoints)**2 + .01**2*np.array(bvrthru)**2) 

plt.errorbar(multqewvls, np.array(multqe)*100, xerr = widths, yerr=multqeErrs*100, color='blue',fmt='o',markersize=8)
print multqe


####Load and average On-sky 2013 throughput, plot up to where data gets messy
path2013 = '/home/srmeeker/ARCONS-pipeline/examples/Pal2013_throughput/'
files = ['115553/hz21_throughput.npz','120100/hz21_throughput.npz','121634/hz21_throughput.npz','122152/hz21_throughput.npz']
fnum=0
allthruput = []
for filename in files:
    t = np.load(path2013+filename)
    wvls2013 = t['wvls']
    thruput = t['throughput']
    if fnum == 0:
        curve2013=thruput
    else:
        curve2013+=thruput
    allthruput.append(thruput)
    fnum+=1
curve2013/=len(files) #get average of each throughput file at each wavelength
#get standard deviation of each file at each wavelength
allthruput = np.array(allthruput)
print allthruput
thruputErrs = np.zeros((len(wvls2013)),dtype=float)
print thruputErrs
for i in xrange(len(wvls2013)):
    thruputErrs[i] = np.std(allthruput[:,i])
print "2013 Thruput errs = ", thruputErrs
#plot 2013 throughput
ax.plot(wvls2013[:-4],curve2013[:-4]*100,linestyle="-.",color='black')#, label=On-sky Total Throughput (2012)')
plt.fill_between(wvls2013[:-4],(curve2013[:-4]-np.abs(thruputErrs[:-4]))*100,(curve2013[:-4]+np.abs(thruputErrs[:-4]))*100,color='PaleGreen',alpha=0.2)

#automatic legend generation
#leg=plt.legend(['Telescope BVR Throughput (2013)','ARCONS QE (2013)', 'On-sky Total Throughput (2013)', 'On-sky Total Throughput (2012)'],'upper right', numpoints=1)
leg=plt.legend(['Telescope BVR\nThroughput','ARCONS QE', "(Telescope TP) *\n(ARCONS QE)", 'Measured Total\nThroughput'],'upper right', numpoints=1,prop={'size':8})

#generate legend manually with right justified text
#vp = leg._legend_box._children[-1]._children[0] 
#for c in vp._children: 
#   c._children.reverse() 
#vp.align="right"

#ax.set_ylim(4E-3,0.04)
ax.set_ylim(4e-3*100,60)
ax.set_xlim(4000,10500)
plt.xlabel(ur"Wavelength (\r{A})")
plt.ylabel(ur"Throughput (\%)")
#set y axis to log scale
ax.set_yscale('log')

# force y-axis formatting to non-sci notation
yax = plt.gca().yaxis 
yax.set_major_formatter(ScalarFormatter()) 
#plt.ticklabel_format(style='plain',useOffset=False)

plt.savefig("FluxCal_LightAccounting.eps",format='eps')

#MAKE THROUGHPUT CURVE FROM absSpectrumG158.py
fig = plt.figure()
ax = fig.add_subplot(111)
#Plot 2012 measured throughput data inverse: flux cal correction curve
ax.plot(wvls,1/(curve/y),'black',linewidth=3)#, label='On-sky Total Throughput (2012)')
ax.set_xlim(4000,11000)
ax.set_ylim(0,250)
plt.xlabel(ur"Wavelength (\r{A})")
plt.ylabel(ur"Spectral Calibration Curve")
plt.savefig("FluxCal_SensitivityCurve.eps",format='eps')

#plt.show()

#np.savez('%s_throughput.npz'%(objectName.strip()),throughput=curve/y,wvls=wvls)

