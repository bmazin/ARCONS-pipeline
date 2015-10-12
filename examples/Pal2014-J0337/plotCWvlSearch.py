import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter,FixedLocator
from util.popup import plotArray,PopUp
from combinePulseProfiles import nSigma
from kuiper.htest import h_fpp
#import figureHeader


#metrics=metrics,wvlLimits=wvlLimits,metricImg=metricImg,fppImg=fppImg,sigmaImg=sigmaImg,wvlPairs=wvlPairs,energyBinWidth=energyBinWidth)
path = '/Scratch/dataProcessing/J0337/realDataCuts.txt'

data = np.loadtxt(path)

wvlStartList = data[:,0]
wvlStopList = data[:,1]
hMetrics = data[:,2]

wvlPairs = zip(wvlStartList,wvlStopList)
wvlLimits = np.unique(np.append(wvlStartList,wvlStopList))
nWvlPairs = len(wvlPairs)
nWvlStarts = len(wvlLimits) - 1
nWvlLimits = len(wvlLimits)

metricImg = np.zeros((nWvlLimits,nWvlLimits))
fppImg = np.zeros((nWvlLimits,nWvlLimits))
sigmaImg = np.zeros((nWvlLimits,nWvlLimits))

iWvlPair = 0
for iWvlStart,wvlStart in enumerate(wvlLimits[:-1]):
    for iWvlEnd,wvlEnd in enumerate(wvlLimits[iWvlStart+1:]):
        assert wvlStart == wvlStartList[iWvlPair]
        assert wvlEnd == wvlStopList[iWvlPair]

        h = hMetrics[iWvlPair]
        pval = h_fpp(h)
        sig = nSigma(1-pval)
        
        metricImg[iWvlStart,iWvlEnd+iWvlStart+1] = h
        fppImg[iWvlStart,iWvlEnd+iWvlStart+1] = pval
        sigmaImg[iWvlStart,iWvlEnd+iWvlStart+1] = sig
        iWvlPair += 1

bestSig = np.max(sigmaImg)
maxH = np.max(metricImg)
maxIndices = np.unravel_index(np.argmax(metricImg),metricImg.shape)
print 'best h at ',maxIndices,' ',wvlLimits[maxIndices[0]],wvlLimits[maxIndices[1]]
print 'h: ',metricImg[maxIndices]
print 'fpp: ',fppImg[maxIndices]
print sigmaImg[maxIndices],'sigmas'


x0 = (23.*11000)/(11000-3600)
c = 3600*x0

def idxToWvl(idx):
    return c/(x0-idx)

def wvlToIdx(wvl):
    return x0-1.*c/wvl

def wvlFormatter(x,p):
    try:
        label = '{:.0f}'.format(idxToWvl(x))
        return label
    except:
        return '0'

yWvlTickLocations = np.append(3600,np.arange(4000,12000,1000))
xWvlTickLocations = np.array([3600,4000,5000,7000,11000])

yIdxTickLocations = np.array([wvlToIdx(wvl) for wvl in yWvlTickLocations])
xIdxTickLocations = np.array([wvlToIdx(wvl) for wvl in xWvlTickLocations])

tickFormatter = FuncFormatter(wvlFormatter)
xTickLocator = FixedLocator(xIdxTickLocations)
yTickLocator = FixedLocator(yIdxTickLocations)

np.savez('metricImg.npz',metricImg=metricImg,wvlLimits=wvlLimits,wvlStartList=wvlStartList,wvlStopList=wvlStopList)

pop = PopUp(showMe=False)
pop.plotArray(metricImg,cbarLabel='h metric',cbarShrink=.82,cbarAspect=15,cbarPad=.05,vmax=13,cmap='hot')
#pop.axes.set_yticks(np.arange(np.shape(wvlLimits)[0]))
pop.axes.tick_params(axis='both', which='major', labelsize=11,labelbottom=True,labeltop=False)
pop.axes.xaxis.set_major_formatter(tickFormatter)
pop.axes.yaxis.set_major_formatter(tickFormatter)
pop.axes.xaxis.set_major_locator(xTickLocator)
pop.axes.yaxis.set_major_locator(yTickLocator)
pop.axes.set_xlabel('upper wavelength limit ($\AA$)')
pop.axes.set_ylabel('lower wavelength limit ($\AA$)')
pop.fig.subplots_adjust(left=.18,right=.93,top=1.0,bottom=None)
pop.fig.set_size_inches(5,4)
#pop.fig.set_tight_layout({})
#pop.fig.tight_layout()
#pop.axes.set_title('H Metric')

#for tick in pop.axes.xaxis.get_major_ticks():
#    tick.label.set_fontsize(11) 
#for tick in pop.axes.yaxis.get_major_ticks():
#    print tick
#    tick.label.set_fontsize(11) 
#    tick.label.set_rotation('vertical')
pop.fig.savefig('metricImg.eps')
pop.show()




