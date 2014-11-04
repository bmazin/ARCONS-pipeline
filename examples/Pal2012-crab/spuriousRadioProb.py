#erfProb.py
#Author: Matt Strader
# this program compares the probabilities of spurious radio giant crab pulses
# calculated  with the assumptions of gaussian white noise vs rfi
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.interpolate
import scipy.stats

phaseShift = 1.-0.681396484375
def indexToPhase(indices):
    radioIndexOffset = 71
    radioArrivalPhases = (indices+71)/2048.+phaseShift
    return radioArrivalPhases
def probsOfGRP(nRadioBins=15,lowerStrengthHistEdge=.12):
    firstDetectedRadioPulseNumber = 3829760547
    lastDetectedRadioPulseNumber = 3829981734

    nRadioPeriodsObserved = lastDetectedRadioPulseNumber-firstDetectedRadioPulseNumber

    labels='iBTOA   pulseNumber BTOA    Noise_Offset    Noise_RMS   Max  Mean Index    TestMax  TestMean TestIndex'
    labels = np.array(labels.split())
    path = '/Scratch/dataProcessing/crabData2/'
    table = np.loadtxt(path+'giantPulseList_P2_3sigma_indices.txt',skiprows=1,usecols=range(len(labels)))
    print 'table',np.shape(table)


    giantDict = dict()
    for iLabel,label in enumerate(labels):
        giantDict[label] = table[:,np.argmax(label==labels)]

    peakToCheck = 2
    noiseStrengthLabel = 'TestMax'
    strengthLabel = 'Max'


    #count number of detections in the test range and the P2 peak range
    #peakDetectionMask = giantDict[strengthLabel]!=0
    peakDetectionMask = giantDict[strengthLabel]>=lowerStrengthHistEdge
    peakDetectionMask2 = giantDict[strengthLabel]>=.13
    print 'peakDetectionMask',np.sum(peakDetectionMask),np.shape(peakDetectionMask)
    print 'peakDetectionMask2',np.sum(peakDetectionMask2),np.shape(peakDetectionMask2)
    noiseDetectionMask = giantDict[noiseStrengthLabel]!=0

    noiseDetections = giantDict[noiseStrengthLabel][noiseDetectionMask]
    peakDetections = giantDict[strengthLabel][peakDetectionMask]
    radioStrength = giantDict[strengthLabel][peakDetectionMask]
    peakDetectedIndices = giantDict['Index'][peakDetectionMask]
    noiseDetectedIndices = giantDict['TestIndex'][noiseDetectionMask]

    clipPeakIndexMask = np.ones(len(peakDetections))
    clipNoiseIndexMask = np.ones(len(noiseDetections))

    startPeakIndex = 1369#np.min(peakDetectedIndices)#1308
    endPeakIndex = 1394#np.max(peakDetectedIndices)#1333
    print 'start,end',startPeakIndex,endPeakIndex
    startNoiseIndex = startPeakIndex-np.min(peakDetectedIndices)+np.min(noiseDetectedIndices)
    endNoiseIndex = endPeakIndex-np.max(peakDetectedIndices)+np.max(noiseDetectedIndices)
    print 'indices min max',np.min(peakDetectedIndices),np.max(peakDetectedIndices)
    clipPeakIndexMask = peakDetectedIndices <= endPeakIndex
    clipPeakIndexMask = np.logical_and(clipPeakIndexMask,peakDetectedIndices >= startPeakIndex)

    clipNoiseIndexMask = noiseDetectedIndices <= endNoiseIndex
    clipNoiseIndexMask = np.logical_and(clipNoiseIndexMask,noiseDetectedIndices >= startNoiseIndex)

    noiseDetections = noiseDetections[clipNoiseIndexMask]
    noiseDetectedIndices = noiseDetectedIndices[clipNoiseIndexMask]
    peakDetections = peakDetections[clipPeakIndexMask]
    peakDetectedIndices = peakDetectedIndices[clipPeakIndexMask]

    nRadioIndexBins = np.max(peakDetectedIndices)-np.min(peakDetectedIndices)+1

    peakDetectionMask3 = giantDict[strengthLabel]>= .175
    noiseDetectionMask3 = giantDict['TestMax']>= .175
    radioStrength3 = giantDict[strengthLabel][peakDetectionMask3]
    peakDetectedIndices3 = giantDict['Index'][peakDetectionMask3]
    noiseDetectedIndices3 = giantDict['TestIndex'][noiseDetectionMask3]
    clipPeakIndexMask3 = peakDetectedIndices3 <= endPeakIndex
    clipPeakIndexMask3 = np.logical_and(clipPeakIndexMask3,peakDetectedIndices3 >= startPeakIndex)
    clipNoiseIndexMask3 = noiseDetectedIndices3 <= endNoiseIndex
    clipNoiseIndexMask3 = np.logical_and(clipNoiseIndexMask3,noiseDetectedIndices3 >= startNoiseIndex)
    peakDetectedIndices3 = giantDict['Index'][peakDetectionMask3]
    noiseDetectedIndices3 = noiseDetectedIndices3[clipNoiseIndexMask3]
    peakDetectedIndices3 = peakDetectedIndices3[clipPeakIndexMask3]

    print 'peakDetectionMask3',np.sum(peakDetectionMask3),np.shape(peakDetectionMask3)
    print 'noiseDetectionMask3',np.sum(noiseDetectionMask3),np.shape(noiseDetectionMask3)
    print 'clipPeakIndexMask3',np.sum(clipPeakIndexMask3),np.shape(clipPeakIndexMask3)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    histPeakIndices,peakIndexBins = np.histogram(peakDetectedIndices3,bins=nRadioIndexBins)
    histNoiseIndices,noiseIndexBins = np.histogram(noiseDetectedIndices3,bins=nRadioIndexBins)
    ax.step(indexToPhase(peakIndexBins[0:-1]),histPeakIndices)
    ax.step(indexToPhase(peakIndexBins[0:-1]),histNoiseIndices)
    print 'nPeak,nNoise',np.sum(histPeakIndices),np.sum(histNoiseIndices)
    ax.set_ylabel('Number of detections')
    ax.set_xlabel('radio phase index')
    ax.set_title('Detections in P2 and TestP2 ranges')


    nDistBins = 70
    peakDist,distBinEdges = np.histogram(peakDetections,bins=nDistBins,range=(lowerStrengthHistEdge,1))
    print 'distBinEdges',distBinEdges
    noiseDist,distBinEdges = np.histogram(noiseDetections,bins=distBinEdges)

#    #plot the histograms in linear scale
#    fig = plt.figure()
#    ax = fig.add_subplot(111)
#    ax.plot(distBinEdges[0:-1],noiseDist,'g',label='noise detections')
#    ax.plot(distBinEdges[0:-1],peakDist-noiseDist,'b',label='\'GRP detections\'')
#    ax.plot(distBinEdges[0:-1],peakDist,'k',label='GRP+noise detections')
#    ax.legend(loc='best')
#    ax.set_xlabel('radio strength')
#    ax.set_title('Detection strength histograms')

    distBinCenters = distBinEdges[0:-1]+np.diff(distBinEdges)/2.

    #now plot histograms in loglog plot.  If the GRP detections have a power law form, it will be a straight line.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog(distBinCenters,noiseDist,'g',label='noise detections')
    ax.loglog(distBinCenters,peakDist-noiseDist,'b',label='\'GRP detections\'')
    ax.loglog(distBinCenters,peakDist,'k',label='GRP+noise detections')
    ax.legend(loc='best')
    ax.set_xlabel('radio strength')
    ax.set_title('Detection strength histograms')

#    print 'log10(strengthBinEdges) log10(nTestP2Detections) log10(nP2Detections) log10(nGRPDetections)=log10(nP2Detections-nTestP2Detections)'
#    for iBin,bin in enumerate(distBinEdges[0:-1]):
#        print iBin,np.log10(bin),np.log10(noiseDist[iBin]),np.log10(peakDist[iBin]),np.log10(peakDist[iBin]-noiseDist[iBin])




    #the probability that a radio detection at a certain radio strength is spurious is the ratio of total false detections in the test range to the total detections (true+false) in the P2 range
    spuriousProb = noiseDist/(1.*peakDist)
    grpProb = 1-spuriousProb
    grpProbF = scipy.interpolate.interp1d(distBinCenters,grpProb,kind='cubic')
    print 'interpolation range',distBinCenters[0],distBinCenters[-1]
    #grpProbAtIndices = grpProb[indices]
    nPoints = 100
    filledStrengthBins = np.linspace(lowerStrengthHistEdge+.01,distBinCenters[-1],nPoints)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(distBinCenters,grpProb,'b',label='P(detection is GRP)')
    ax.plot(distBinCenters,spuriousProb,'g',label='P(detection is noise)')
    ax.set_xlabel('giant radio strength')
    ax.legend(loc='upper right')

    
    #only include GRP that fall in the desired index range and strength range (threshold)

    radioStrength2 = giantDict[strengthLabel][peakDetectionMask2]
    peakDetectedIndices2 = giantDict['Index'][peakDetectionMask2]
    clipPeakIndexMask2 = peakDetectedIndices2 <= endPeakIndex
    clipPeakIndexMask2 = np.logical_and(clipPeakIndexMask2,peakDetectedIndices2 >= startPeakIndex)

    #percentiles = np.arange(0,100.1,100./nRadioBins)
    #percentiles = np.logspace(-1,2,nRadioBins,base=2)
    #percentiles = np.array([0,15,30,42,47,52,57,62,67,75,85,90])
    percentiles = np.array([0,15,30,44,50,55,60,65,70,75,85,90])
    #percentiles = np.array([0,18,37,53,62,70,77,84,90])
    #percentiles = np.array([0,15,34,50,59,66,73,80,86,90])
    #percentiles = np.array([0,15,30,45,50,55,60,65,70,75,85,90])
    radioStrengthBins = []
    for perc in percentiles:
        radioStrengthBins.append(scipy.stats.scoreatpercentile(radioStrength2[clipPeakIndexMask2],perc))

    #radioStrengthBins = np.linspace(.13,.3,nRadioBins)

    radioStrengthBinCenters = radioStrengthBins[0:-1]+np.diff(radioStrengthBins)/2.
    #grpProbAtIndices = grpProbF(radioStrengthBinCenters)
    #ax.plot(radioStrengthBins[0:len(indices)],grpProbAtIndices,'r.')
    ax.set_xlim([.12,.3])

#    radioStrengthCutoff = .175
#    radioCutoffMask = radioStrength >= radioStrengthCutoff
#    #only include GRP that fall in the desired index range and strength range (threshold)
#    radioPeakMask = np.logical_and(clipPeakIndexMask,radioCutoffMask)

    #a different sized mask is needed for strength>=.13

    return {'radioPeakMask':clipPeakIndexMask2,'radioBins':radioStrengthBins,'grpProbFunc':grpProbF,'grpProbFuncRange':[distBinCenters[0],distBinCenters[-1]]}




if __name__=='__main__':
    binsString = """0.130001     0.13137848   0.13254296   0.13372344   0.13481184
       0.1358914    0.13707388   0.13831236   0.13966668   0.14120796
          0.1428228    0.14459556   0.14662752   0.14906848   0.15211916
             0.1555064    0.160049     0.16608192   0.17604684   0.19213508
                0.222151     0.26999828   0.3452324    0.477951     0.75829176
                  23.311523 """
    radioStrengthBins = np.array([np.double(s) for s in binsString.split()])
    probsOfGRP()
    plt.show()
