import numpy as np
import os
import matplotlib.pyplot as plt
from util.popup import PopUp,plotArray,pop
import itertools

from fitFunctions import gaussian
import mpfit
from peakWidth import peakWidth
import pickle
from astropy.stats.funcs import sigma_clip


#from http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

#from http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


def fitGauss(dataList,nBins=201):

    hist,histBinEdges = np.histogram(dataList,bins=nBins)
    histBinCenters = histBinEdges[0:-1]+np.diff(histBinEdges)/2.

    amplitude = 0.95*np.max(hist)
    x_offset = histBinCenters[np.argmax(hist)]
    sigma = np.std(dataList)*.8
    y_offset = 1.e-8

    params=[sigma, x_offset, amplitude, y_offset]  # First guess at fit params

    #errs = np.sqrt(hist) #for poisson distributed data
    #errs[np.where(errs == 0.)] = 1.

    nListVals = 1.*np.sum(hist)
    pvals = hist/nListVals
    #assume errors are the multinomial distribution
    errs = np.sqrt(nListVals*pvals*(1-pvals))
    errs[np.where(errs == 0.)] = 1.


    quiet = True

    parinfo = [ {'n':0,'value':params[0],'limits':[sigma/10., 10*sigma], 'limited':[True,True],'fixed':False,'parname':"Sigma",'error':0},
       {'n':1,'value':params[1],'limits':[x_offset-sigma*2, x_offset+sigma*2],'limited':[True,True],'fixed':False,'parname':"x offset",'error':0},
       {'n':2,'value':params[2],'limits':[0.5*amplitude, 2.*amplitude],'limited':[True,True],'fixed':False,'parname':"Amplitude",'error':0},
       {'n':3,'value':params[3],'limited':[False,False],'fixed':True,'parname':"y_offset",'error':0}]

    fa = {'x':histBinCenters,'y':hist,'err':errs}

    m = mpfit.mpfit(gaussian, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
    if m.status <= 0:
        print m.status, m.errmsg

    mpp = m.params                                #The fit params
    mpperr = m.perror

    for k,p in enumerate(mpp):
        parinfo[k]['value'] = p
        parinfo[k]['error'] = mpperr[k]
        if k==0: sigma = p
        if k==1: x_offset = p
        if k==2: amplitude = p
        if k==3: y_offset = p

    def gaussFitFunc(x):
        return y_offset + amplitude * np.exp( - (( x - x_offset)**2) / ( 2. * (sigma**2)))
    gaussFit = gaussFitFunc(histBinCenters)
    resolution = np.abs(x_offset/(2.355*sigma))

    return {'gaussFit':gaussFit,'resolution':resolution,'sigma':sigma,'x_offset':x_offset,'amplitude':amplitude,'y_offset':y_offset,'hist':hist,'histBinEdges':histBinEdges,'gaussFitFunc':gaussFitFunc,'histBinCenters':histBinCenters,'parinfo':parinfo}

flatTstamp = '20121211'
flatType = 'flat'
#flatTstamp = '20121212-074700' #Geminga sky flat
folder = '/Scratch/dataProcessing/flatTests/'
lowerWvlCut = 4000
upperWvlCut = 8000
fileName = '{}Cube_{}m_wvl{}-{}.npz'.format(flatType,flatTstamp,lowerWvlCut,upperWvlCut)
#fileName = 'flatCube_{}.npz'.format(flatTstamp)
print fileName

imgDict = np.load(os.path.join(folder,fileName))
beforeImg = imgDict['beforeImg']
afterImg = imgDict['afterImg']

bStableMask = False
if bStableMask:
    sdevImg = np.load('stableMask.npz')['sdevImg']
    sdevThreshold = 4.
    stableMask = np.logical_and(100.*sdevImg<sdevThreshold,sdevImg>0.)
    beforeImg[~stableMask] = 0
    afterImg[~stableMask] = 0
    print 'masked stable sdev<{}%'.format(sdevThreshold)

bCoolMask = False
if bCoolMask:
    effIntImg = np.load('stableMask.npz')['effIntImg']
    effIntThreshold = .999
    mask = effIntImg >= effIntThreshold
    beforeImg[~mask] = 0
    afterImg[~mask] = 0
    print 'masked effIntTime fraction >{}'.format(effIntThreshold)

nRows,nCols = np.shape(beforeImg)

clippedBeforeImg = sigma_clip(beforeImg,sig=3)
clipMask = clippedBeforeImg.mask
beforeImg[clipMask] = 0
afterImg[clipMask] = 0

deadBeforeImg = (beforeImg == 0)
deadAfterImg = (afterImg == 0)

beforeList = beforeImg[beforeImg != 0]
afterList = afterImg[afterImg != 0]

afterImg[deadAfterImg] = np.nan
beforeImg[deadBeforeImg] = np.nan

plotArray(title='without flatcal',image=beforeImg)
#plotArray(title='with flatcal',image=afterImg)

def plotFunc(fig,axes):
    axes.plot(beforeHistEdges[0:-1],beforeHist,label='without flatcal')
    axes.plot(afterHistEdges[0:-1],afterHist,label='with flatcal')
    axes.set_title('Distribution of pixel counts')
    axes.set_xlabel('Counts')
    axes.set_ylabel('Num of Pixels')
    axes.legend()
#pop(plotFunc=plotFunc)

print 'before count',len(beforeList)
print 'after count',len(afterList)

# Fit a 3rd order, 2d polynomial to the non-flatcal image
xx,yy = np.meshgrid(np.arange(nCols),np.arange(nRows))

z = beforeImg.ravel()
x = xx.ravel()
y = yy.ravel()

x = x[z != 0]
y = y[z != 0]
z = z[z != 0]

beforePolyFitCoeffs = polyfit2d(x,y,z)

# Evaluate it on a grid...
beforeImgFit = polyval2d(xx, yy, beforePolyFitCoeffs)


# Fit a 3rd order, 2d polynomial to the flatcal image

z = afterImg.ravel()
x = xx.ravel()
y = yy.ravel()

x = x[z != 0]
y = y[z != 0]
z = z[z != 0]

afterPolyFitCoeffs = polyfit2d(x,y,z)

# Evaluate it on a grid...
afterImgFit = polyval2d(xx, yy, afterPolyFitCoeffs)

plotArray(beforeImgFit,vmin=0,title='poly fit to non-flatcal\'d image')
plotArray(afterImgFit,vmin=0,title='poly fit to flatcal\'d image')

#afterImgSub = afterImg / np.mean(afterList)
afterImgSub = np.array(afterImg)
#afterImgSub[deadAfterImg] = 0

beforeImgSub = beforeImg * np.mean(beforeImgFit) / beforeImgFit
#beforeImgSub[deadBeforeImg] = 0

plotArray(beforeImgSub,title='unflatcal\'d, scaled by gray illumination')
plotArray(afterImgSub,title='flatcal\'d')


subAfterList = afterImgSub[~np.isnan(afterImgSub)]
subBeforeList = beforeImgSub[~np.isnan(beforeImgSub)]
nBins = int(len(afterList)/6.)
subBeforeHist,subBeforeHistEdges = np.histogram(subBeforeList,bins=nBins)
subAfterHist,subAfterHistEdges = np.histogram(subAfterList,bins=subBeforeHistEdges)
beforeHist,beforeHistEdges = np.histogram(beforeList,bins=nBins)
afterHist,afterHistEdges = np.histogram(afterList,bins=beforeHistEdges)

rawFwhm = peakWidth(beforeHistEdges[0:-1],beforeHist)
beforeFwhm = peakWidth(subBeforeHistEdges[0:-1],subBeforeHist)
afterFwhm = peakWidth(subAfterHistEdges[0:-1],subAfterHist)

print 'raw sigma', rawFwhm['sigma'],'{:.1f}%'.format(100.*rawFwhm['sigma']/rawFwhm['peakX']),np.sum((beforeList < 3200) | (beforeList > 4000))
print 'gray scaled sigma', beforeFwhm['sigma'],'{:.1f}%'.format(100.*beforeFwhm['sigma']/beforeFwhm['peakX']),np.sum((subBeforeList < 3200) | (subBeforeList > 4000))
print 'flatfielded sigma', afterFwhm['sigma'],'{:.1f}%'.format(100.*afterFwhm['sigma']/afterFwhm['peakX']),np.sum((subAfterList < 3200) | (subAfterList > 4000))

shotNoiseSigma = np.sqrt(afterFwhm['peakX'])
print 'shot noise sigma',shotNoiseSigma,'{:.1f}%'.format(100.*shotNoiseSigma/afterFwhm['peakX'])

print 'noise ratio',afterFwhm['sigma']/shotNoiseSigma

# Plot

def plotFunc(fig,axes):
    axes.set_title('Distribution of pixel counts')
    axes.plot(beforeHistEdges[0:-1],beforeHist,color='k',label='raw')
    axes.plot(subBeforeHistEdges[0:-1],subBeforeHist,color='b',label='grayscale corrected')
    axes.plot(subAfterHistEdges[0:-1],subAfterHist,'r',label='flatfield applied')
    axes.set_xlabel('Counts')
    axes.set_ylabel('Num of Pixels')
    axes.legend()
pop(plotFunc=plotFunc)

plotDict = {'rawHistEdges':beforeHistEdges,'rawHist':beforeHist,'illumHistEdges':subBeforeHistEdges,'illumHist':subBeforeHist,'flatHistEdges':subAfterHistEdges,'flatHist':subAfterHist}
pickle.dump(plotDict,open( 'flatHist.pickle', 'wb' ))
print plotDict['rawHist']
print plotDict['illumHist']
print plotDict['flatHist']

#subBeforeGaussFitDict = fitGauss(subBeforeList)
#histBinEdges = subBeforeGaussFitDict['histBinEdges']
#histBinCenters = subBeforeGaussFitDict['histBinCenters']
#x = np.linspace(np.min(histBinEdges),np.max(histBinEdges),1000)
#y = subBeforeGaussFitDict['gaussFitFunc'](x)
#
#subAfterGaussFitDict = fitGauss(subAfterList)
#histBinEdges = subAfterGaussFitDict['histBinEdges']
#histBinCenters = subAfterGaussFitDict['histBinCenters']
#
#x2 = np.linspace(np.min(histBinEdges),np.max(histBinEdges),1000)
#y2 = subAfterGaussFitDict['gaussFitFunc'](x2)
#
#def plotFunc(fig,axes):
#    hist = subBeforeGaussFitDict['hist']
#    axes.plot(x,y,'b')
#    axes.errorbar(subBeforeGaussFitDict['histBinCenters'],hist,color='m')
#    axes.plot(x2,y2,'g')
#    hist = subAfterGaussFitDict['hist']
#    axes.errorbar(subAfterGaussFitDict['histBinCenters'],hist,color='r')
#    axes.set_title('Distribution of pixel counts')
#pop(plotFunc=plotFunc)
#
#beforeSigmaDict = (item for item in subBeforeGaussFitDict['parinfo'] if item['parname'] == 'Sigma').next()
#print beforeSigmaDict
#print 'non-flatcal\'d sigma',beforeSigmaDict['value'],'+/-',beforeSigmaDict['error']
#afterSigmaDict = (item for item in subAfterGaussFitDict['parinfo'] if item['parname'] == 'Sigma').next()
#print 'flatcal\'d sigma',afterSigmaDict['value'],'+/-',afterSigmaDict['error']

