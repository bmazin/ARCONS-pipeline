import numpy as np
import matplotlib.pyplot as plt
from util import utils
from util.readDict import readDict
from scipy import pi
from mpltools	import style
import figureHeader

# try new plotting style for pipeline paper
#style.use('ggplot')

####LOAD COROT18 LIGHTCURVE######
param = readDict()
param.read_from_file('corot18bparams.dict')

FileName = param['npzLoadFitFile']
FramesPerFile = param['FramesPerFile']
print FramesPerFile
TotalNumFiles = param['TotalNumFiles']
NumFrames = FramesPerFile*TotalNumFiles
IntTime = param['integrationTime']

#FoldPeriod = 0.01966127 #This is in fractions of a day
t = np.load(FileName)
params = t['params']
jd = t['jd']
#params[:,0] are the height offset of the gaussian (so the number of background pulses)
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
#jd2 = (jd/FoldPeriod)%1.

#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(jd, xpos,'k.')
#ax.plot(jd, ypos,'k.',color='red')
#plt.show()

Corotcurve = 2*pi*amps*widths**2
#need to convert to Jy: curve*wvlrange*aveE/(1Hz*intTime*Area)

fwhm = 2*np.sqrt(2*np.log(2))*widths#pixels
#fwhm = 0.5*fwhm #arcsec
#fwhm = widths
corotWidths = fwhm

CorotmeanXpos = utils.mean_filterNaN(xpos,size=7)
CorotmeanYpos = utils.mean_filterNaN(ypos,size=7)
#jd = jd[curve<=500]
#curve = curve[curve<=500]
#jd = jd[curve>=100]
#curve = curve[curve>=100]


####LOAD COMPANION LIGHT CURVE######
param = readDict()
param.read_from_file('corot18Compparams.dict')

FileName = param['npzLoadFitFile']
FramesPerFile = param['FramesPerFile']
print FramesPerFile
TotalNumFiles = param['TotalNumFiles']
NumFrames = FramesPerFile*TotalNumFiles
IntTime = param['integrationTime']

#FoldPeriod = 0.01966127 #This is in fractions of a day
t = np.load(FileName)
params = t['params']
jd = t['jd']
#params[:,0] are the height offset of the gaussian (so the number of background pulses)
amps = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]
#jd2 = (jd/FoldPeriod)%1.

#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(jd, xpos,'k.')
#ax.plot(jd, ypos,'k.',color='red')
#plt.show()


#fig = plt.figure()
#ax = fig.add_subplot(111)

Compcurve = 2*pi*amps*widths**2
#need to convert to Jy: curve*wvlrange*aveE/(1Hz*intTime*Area)

fwhm = 2*np.sqrt(2*np.log(2))*widths#pixels
#fwhm = 0.5*fwhm #arcsec
#fwhm = widths
compWidths = fwhm

CompmeanXpos = utils.mean_filterNaN(xpos,size=7)
CompmeanYpos = utils.mean_filterNaN(ypos,size=7)
#jd = jd[curve<=500]
#curve = curve[curve<=500]
#jd = jd[curve>=100]
#curve = curve[curve>=100]


#ax.plot(jd, corotWidths,'k.',color='black')
#ax.plot(jd, compWidths,'k.',color='red')
#plt.show()


fig = plt.figure()
ax = fig.add_subplot(111)

#curve/=np.median(curve)
#fwhm/=np.median(fwhm)
ax.plot(jd[jd<np.median(jd)]-2456269,Corotcurve[jd<np.median(jd)],'k.',color='black',markersize=8)
#ax.plot(jd[jd<np.median(jd)]-2456269,Compcurve[jd<np.median(jd)],'k.',color='red')
ax.set_xlabel("JD-2456269")#,fontsize=20)
ax.set_ylabel("Counts")#,fontsize=20)
#ax.patch.set_visible(False)
#ax.set_title("Uncorrected Corot-18 Lightcurve")
#ax.tick_params(axis='both', which='major', labelsize=14)
#plt.ylim([80000/180000,1])
#plt.ylim([80000/3.0,180000/3.0])
#plt.yscale('log')
plt.savefig("Corot_all.eps",format='eps',bbox_inches='tight')
#plt.show()

subCounts = Corotcurve[(jd>2456269.8077) & (jd<2456269.8141)]
subJD = jd[(jd>2456269.8077) & (jd<2456269.8141)]
meanCounts = np.ones(len(subJD),dtype=float)
meanCounts *= np.mean(subCounts)

meanVal = np.mean(subCounts)
print "mean = %i"%meanVal
meanCounts -=meanVal
subCounts -=meanVal
stdev = np.std(subCounts)
print "stdev = %f"%stdev
print "percent = %f %%"%(stdev/meanVal *100)

fig = plt.figure()
ax2 = fig.add_subplot(111)

#ax2 = plt.axes([.4,.15,.5,.5])
#ax2.set_xlabel("JD-2456269",fontsize=20)
ax2.set_ylabel("Counts-mean",fontsize=20)
ax2.tick_params(axis='both', which='major', labelsize=12)
ax2.plot(subJD-2456269,subCounts,'k.',markersize=12)
ax2.plot(subJD-2456269,meanCounts,color='grey',linewidth=2)
ax2.axhspan(-stdev,stdev,facecolor='PowderBlue',alpha=0.2)
ax2.set_xlim([.808,.814])
plt.savefig("Corot_zoom.eps",format='eps')
#plt.show()


ratio = Corotcurve/Compcurve
jd = jd[ratio!=np.inf]
ratio = ratio[ratio!=np.inf]
jd = jd[ratio!=np.nan]
ratio = ratio[ratio!=np.nan]
jd = jd[ratio<=1000]
ratio = ratio[ratio<=1000]

#fig = plt.figure()
#ax = fig.add_subplot(111)
#ax.plot(jd,ratio,'k.')
#plt.show()

