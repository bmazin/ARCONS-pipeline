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

if __name__=='__main__':
    flatLabel = '20121211'
    #flatLabel = '20121212-074700'
    lowerWvlCut = 4000
    upperWvlCut = 11000
    path = '/Scratch/dataProcessing/flatTests/'
    outfilename = os.path.join(path,'{}_{}-{}.txt'.format(flatLabel,lowerWvlCut,upperWvlCut))
    #outfile = open(outfilename,'a')

    #date = '20121210'
    #obsTimestamp = '20121211-051650' #hr9087
    #obsTimestamp = '20121211-124809' #J0926
    #obsTimestamp = '20121211-125312' #J0926
    #obsTimestamp = '20121211-125814' #J0926
    #obsTimestamp = '20121211-130316' #J0926
    #obsTimestamp = '20121211-130818' #J0926
    #obsTimestamp = '20121211-131320' #J0926
    #obsTimestamp = '20121211-131822' #J0926
    #obsTimestamp = '20121211-132324' #last J0926
    #obsTimestamp = '20121211-134751' # feige 66 sky
    #obsTimestamp = '20121211-134812' # feige 66 sky
    #obsTimestamp = '20121211-134914' # feige 66 sky
    #obsTimestamp = '20121211-135016' #early twilight
    #obsTimestamp = '20121211-135118' #later twilight
    #obsTimestamp = '20121211-135220' #later twilight
    #obsTimestamp = '20121211-135322' #later twilight
    #obsTimestamp = '20121211-135424' #later twilight
    #obsTimestamp = '20121211-135526' #later twilight
    #obsTimestamp = '20121211-135628' #late twilight
    #obsTimestamp = '20121211-135730' #late twilight


    #date = '20121211'
    #obsTimestamp = '20121212-085730' #geminga
    #obsTimestamp = '20121212-090233' #geminga
    #obsTimestamp = '20121212-090735' #geminga
    #obsTimestamp = '20121212-091237' #geminga

    #obsTimestamp = '20121212-095809' #ptfo 8-8695
    #obsTimestamp = '20121212-100311' #ptfo 8-8695
    #obsTimestamp = '20121212-100813' #ptfo 8-8695
    #obsTimestamp = '20121212-101315' #ptfo 8-8695
    #obsTimestamp = '20121212-101817' #ptfo 8-8695

    #obsTimestamp = '20121212-103030' #geminga
    #obsTimestamp = '20121212-104035' #geminga
    #obsTimestamp = '20121212-105326' #geminga
    #obsTimestamp = '20121212-110330' #geminga
    #obsTimestamp = '20121212-111334' #geminga

    #obsTimestamp = '20121212-112709' #J0926
    #obsTimestamp = '20121212-113212' #J0926
    #obsTimestamp = '20121212-113714' #J0926
    #obsTimestamp = '20121212-114216' #J0926
    #obsTimestamp = '20121212-114718' #J0926
    #obsTimestamp = '20121212-115220' #J0926
    #obsTimestamp = '20121212-115722' #J0926
    #obsTimestamp = '20121212-120224' #J0926
    #obsTimestamp = '20121212-120727' #J0926
    #obsTimestamp = '20121212-121229' #J0926
    #obsTimestamp = '20121212-121732' #J0926
    #obsTimestamp = '20121212-122234' #J0926
    #obsTimestamp = '20121212-122736' #J0926
    #obsTimestamp = '20121212-123238' #J0926
    #obsTimestamp = '20121212-123740' #J0926
    #obsTimestamp = '20121212-124242' #J0926
    #obsTimestamp = '20121212-124744' #J0926
    #obsTimestamp = '20121212-125246' #J0926
    #obsTimestamp = '20121212-125748' #J0926
    #obsTimestamp = '20121212-130250' #J0926
    #obsTimestamp = '20121212-130752' #J0926
    #obsTimestamp = '20121212-131254' #J0926
    #obsTimestamp = '20121212-131756' #J0926
    #obsTimestamp = '20121212-132258' #J0926
    #obsTimestamp = '20121212-132800' #J0926
    #obsTimestamp = '20121212-133303' #J0926
    #obsTimestamp = '20121212-134024' #early twilight
    #obsTimestamp = '20121212-134127' #twilight
    #obsTimestamp = '20121212-134229' #twilight
    #obsTimestamp = '20121212-134331' #twilight
    #obsTimestamp = '20121212-134433' #twilight
    #obsTimestamp = '20121212-134535' #twilight###
    #obsTimestamp = '20121212-134637' #twilight
    #obsTimestamp = '20121212-134739' #twilight
    #obsTimestamp = '20121212-134841' #twilight
    #obsTimestamp = '20121212-134943' #twilight
    #obsTimestamp = '20121212-135045' #twilight
    #obsTimestamp = '20121212-135147' #twilight
    #obsTimestamp = '20121212-135249' #twilight
    #obsTimestamp = '20121212-135351' #twilight
    #obsTimestamp = '20121212-135453' #twilight
    #obsTimestamp = '20121212-135555' #twilight
    #obsTimestamp = '20121212-135657' #twilight
    #obsTimestamp = '20121212-135759' #twilight
    #obsTimestamp = '20121212-135901' #twilight
    #obsTimestamp = '20121212-140003' #twilight
    #obsTimestamp = '20121212-140105' #late twilight

    #date = '20140923'
    #obsTimestamp = '20140924-065535' #early 1SWASP sky
    #obsTimestamp = '20140924-123738' #early twilight
    #obsTimestamp = '20140924-125220' #late twilight

    #date = '20140924'
    #obsTimestamp ='20140925-030054' #G24-9 sky, faint objects on array
    #obsTimestamp ='20140925-052248' #V407 Vul
    #obsTimestamp ='20140925-052804' #V407 Vul
    #obsTimestamp ='20140925-053307' #V407 Vul
    #obsTimestamp ='20140925-053810' #V407 Vul

    #obsTimestamp ='20140925-110513' #psr J0337
    #obsTimestamp ='20140925-111029' #psr J0337
    #obsTimestamp ='20140925-111532' #psr J0337
    #obsTimestamp ='20140925-112035' #psr J0337
    #obsTimestamp ='20140925-112538' #psr J0337
    #obsTimestamp ='20140925-113041' #psr J0337
    #obsTimestamp ='20140925-113544' #psr J0337
    #obsTimestamp ='20140925-114047' #psr J0337
    #obsTimestamp ='20140925-114550' #psr J0337
    #obsTimestamp ='20140925-115053' #mid psr J0337
    #obsTimestamp ='20140925-115605' #psr J0337
    #obsTimestamp ='20140925-120112' #psr J0337
    #obsTimestamp ='20140925-120615' #psr J0337
    #obsTimestamp ='20140925-121118' #psr J0337
    #obsTimestamp ='20140925-121621' #psr J0337
    #obsTimestamp ='20140925-122124' #psr J0337
    #obsTimestamp ='20140925-122627' #psr J0337
    #obsTimestamp ='20140925-123130' #psr J0337
    #obsTimestamp ='20140925-123633' #last PSR J0337 file for 20140924
    #obsTimestamp ='20140925-124254' #early twilight
    #obsTimestamp ='20140925-124357' #early twilight
    #obsTimestamp ='20140925-124500' #early twilight
    #obsTimestamp ='20140925-124603' #early twilight
    #obsTimestamp ='20140925-124706' #early twilight
    #obsTimestamp ='20140925-124809' #early twilight
    #obsTimestamp ='20140925-124912' #early twilight
    #obsTimestamp ='20140925-125015' #late twilight
    #obsTimestamp ='20140925-125118' #late twilight
    #obsTimestamp ='20140925-125221' #late twilight
    #obsTimestamp ='20140925-125324' #late twilight
    #obsTimestamp ='20140925-125427' #late twilight


    #date = '20140925'
    #obsTimestamp ='20140926-104528' #psr J0337

    #flatTstamp = '20121212-074700' #Geminga sky flat
    folder = '/Scratch/dataProcessing/flatTests/'
    fileName = 'flattenedCubes_obs{}_flat{}_wvl{}-{}.npz'.format(obsTimestamp,flatLabel,lowerWvlCut,upperWvlCut)
    print fileName

    imgDict = np.load(os.path.join(folder,fileName))
    rawImg = imgDict['rawImg']
    flatImg = imgDict['flatImg']
    illumImg = imgDict['illumImg']
    exptime = imgDict['exptime']
    nRows,nCols = np.shape(flatImg)

    nanMask = np.isnan(rawImg) | np.isnan(flatImg) | np.isnan(illumImg)

    rawImg = np.ma.array(rawImg,mask=nanMask)
    flatImg = np.ma.array(flatImg,mask=nanMask)
    illumImg = np.ma.array(illumImg,mask=nanMask)

    bClipData = True
    sig=3.
    if bClipData:
        iters=None
        cenfunc = np.ma.median
        rawImg = sigma_clip(rawImg,sig=sig,iters=iters,cenfunc=cenfunc)
        flatImg = sigma_clip(flatImg,sig=sig,iters=iters,cenfunc=cenfunc)
        illumImg = sigma_clip(illumImg,sig=sig,iters=iters,cenfunc=cenfunc)
        
        clipMask = rawImg.mask | flatImg.mask | illumImg.mask
        rawImg.mask = clipMask
        flatImg.mask = clipMask
        illumImg.mask = clipMask
        
        print 'removed ',np.sum(clipMask)-np.sum(nanMask),'pixels'

    rawList = rawImg[~rawImg.mask]
    flatList = flatImg[~flatImg.mask]
    illumList = illumImg[~illumImg.mask]

    #plotArray(title='with flatcal',image=flatImg)

    print 'raw count',len(rawList)
    print 'flat count',len(flatList)
    print 'illum count',len(illumList)

    plotArray(rawImg,title='raw')
    plotArray(illumImg,title='illumination corrected')
    plotArray(flatImg,title='flatfield corrected')

    binWidth = 20 #counts
    nBins = int(1.*len(flatList)/binWidth)
    rawHist,rawHistEdges = np.histogram(rawList,bins=nBins)
    flatHist,flatHistEdges = np.histogram(flatList,bins=rawHistEdges)
    illumHist,illumHistEdges = np.histogram(illumList,bins=rawHistEdges)

    rawFwhm = peakWidth(rawHistEdges[0:-1],rawHist)
    flatFwhm = peakWidth(flatHistEdges[0:-1],flatHist)
    illumFwhm = peakWidth(illumHistEdges[0:-1],illumHist)

    avgCounts = int(flatFwhm['peakX'])
    rawFwhmPercent = 100.*rawFwhm['sigma']/avgCounts
    illumFwhmPercent = 100.*illumFwhm['sigma']/avgCounts
    flatFwhmPercent = 100.*flatFwhm['sigma']/avgCounts
    shotNoiseSigma = np.sqrt(flatFwhm['peakX'])
    shotPercent = 100.*shotNoiseSigma/avgCounts
    flatToShotNoiseRatio = flatFwhm['sigma']/shotNoiseSigma

    print 'avg counts',avgCounts
    print 'raw sigma', rawFwhm['sigma'],'{:.1f}%'.format(rawFwhmPercent)

    print 'illumination corrected sigma', illumFwhm['sigma'],
    print '{:.1f}%'.format(illumFwhmPercent)

    print 'flatfield corrected sigma', flatFwhm['sigma'],
    print '{:.1f}%'.format(flatFwhmPercent)

    print 'shot noise sigma',shotNoiseSigma,'{:.1f}%'.format(shotPercent)

    print 'noise ratio',flatToShotNoiseRatio

    # Plot

    def plotFunc(fig,axes):
        axes.set_title('Distribution of pixel counts')
        axes.plot(rawHistEdges[0:-1],rawHist,color='k',label='raw')
        axes.plot(illumHistEdges[0:-1],illumHist,color='b',label='illum')
        axes.plot(flatHistEdges[0:-1],flatHist,'r',label='flatfield')
        axes.set_xlabel('Counts')
        axes.set_ylabel('Num of Pixels')
        axes.legend()
    pop(plotFunc=plotFunc)

    plotDict = {'rawHistEdges':rawHistEdges,'rawHist':rawHist,'illumHistEdges':illumHistEdges,'illumHist':illumHist,'flatHistEdges':flatHistEdges,'flatHist':flatHist,'flatLabel':flatLabel,'lowerWvlCut':lowerWvlCut,'upperWvlCut':upperWvlCut,'imgFilename':outfilename,'obsDate':date,'obsTimestamp':obsTimestamp,'bClippedImgs':bClipData,'nBins':nBins}
    pickle.dump(plotDict,open( 'flatHist.pickle', 'wb' ))
    #outfile.write('{}\t{}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.1f}\t{:.2f}\t{}\n'.format(obsTimestamp,avgCounts,rawFwhmPercent,illumFwhmPercent,flatFwhmPercent,shotPercent,flatToShotNoiseRatio,exptime))
    
    


