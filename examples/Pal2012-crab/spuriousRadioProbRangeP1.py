#erfProb.py
#Author: Matt Strader
# this program compares the probabilities of spurious radio giant crab pulses
# calculated  with the assumptions of gaussian white noise vs rfi
import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.interpolate
import scipy.stats

phaseShift = 1.-0.677001953125#found with findOpticalPeak.p
def indexToPhase(indices):
    radioIndexOffset = 0.5
    radioArrivalPhases = (indices+radioIndexOffset)/2048.+phaseShift
    return radioArrivalPhases
def probsOfGRP(nRadioBins=15,strengthMin=.175,startPeakIndex=None,endPeakIndex=None):
    firstDetectedRadioPulseNumber = 3829760547
    lastDetectedRadioPulseNumber = 3829981734

    nRadioPeriodsObserved = lastDetectedRadioPulseNumber-firstDetectedRadioPulseNumber

    labels='iBTOA   pulseNumber BTOA    Noise_Offset    Noise_RMS   Max  Mean Index    TestMax  TestMean TestIndex'
    labels = np.array(labels.split())
    path = '/Scratch/dataProcessing/crabData2/'
    table = np.loadtxt(path+'giantPulseList_P1_3sigma_indices.txt',skiprows=1,usecols=range(len(labels)))
    #labels='Noise_Offset Noise_RMS Max_P1 Mean_P1 Prob_P1 TestMax_P1 TestMean_P1 TestProb_P1 Max_P2 Mean_P2 Prob_P2 TestMax_P2 TestMean_P2 TestProb_P2'
    #labelsOld = 'BTOA Max_P1 Mean_P1 Max_P2 Mean_P2'
    #labelsOld2 = 'BTOA Max_P1 Mean_P1 Max_P2 Mean_P2'

    #mask out any interpulses in pulses that also have a (weak) main pulse GRP
    overlapPNs = np.load('overlapP1.npz')['overlap']
    pulseNumbers = table[:,np.argmax('pulseNumber'==labels)]
    mainPulseMask = np.logical_not(np.in1d(pulseNumbers,overlapPNs))
    #mainPulseMask = np.logical_or(np.logical_not(mainPulseMask),table[:,np.argmax('TestMax'==labels)]>0.)
    table = table[mainPulseMask]
    print 'table',len(table)

    giantDict = dict()
    for iLabel,label in enumerate(labels):
        giantDict[label] = table[:,np.argmax(label==labels)]

    noiseStrengthLabel = 'TestMax'
    strengthLabel = 'Max'


    nDistBins = 50
    #count number of detections in the test range and the P2 peak range
    peakDetectionMask = giantDict[strengthLabel]>=strengthMin
    noiseDetectionMask = giantDict[noiseStrengthLabel]>=strengthMin
    print 'peakDetectionMask',np.sum(peakDetectionMask),np.shape(peakDetectionMask)

    noiseDetections = giantDict[noiseStrengthLabel][noiseDetectionMask]
    peakDetections = giantDict[strengthLabel][peakDetectionMask]
    radioStrength = giantDict[strengthLabel][peakDetectionMask]
    peakDetectedIndices = giantDict['Index'][peakDetectionMask]
    noiseDetectedIndices = giantDict['TestIndex'][noiseDetectionMask]

    clipPeakIndexMask = np.ones(len(peakDetections))
    clipNoiseIndexMask = np.ones(len(noiseDetections))

    if startPeakIndex == None:
        startPeakIndex = np.min(peakDetectedIndices)#1308
    if endPeakIndex == None:
        endPeakIndex = np.max(peakDetectedIndices)#1333
    print 'start,end',startPeakIndex,endPeakIndex,indexToPhase(startPeakIndex),indexToPhase(endPeakIndex)
    startNoiseIndex = startPeakIndex-np.min(peakDetectedIndices)+np.min(noiseDetectedIndices)
    endNoiseIndex = endPeakIndex-np.max(peakDetectedIndices)+np.max(noiseDetectedIndices)
    print 'peak indices min max',np.min(peakDetectedIndices),np.max(peakDetectedIndices)
    print 'noise indices min max',np.min(noiseDetectedIndices),np.max(noiseDetectedIndices)
    clipPeakIndexMask = peakDetectedIndices <= endPeakIndex
    clipPeakIndexMask = np.logical_and(clipPeakIndexMask,peakDetectedIndices >= startPeakIndex)

    clipNoiseIndexMask = noiseDetectedIndices <= endNoiseIndex
    clipNoiseIndexMask = np.logical_and(clipNoiseIndexMask,noiseDetectedIndices >= startNoiseIndex)
    print np.sum(clipPeakIndexMask),np.sum(clipNoiseIndexMask)

    noiseDetections = noiseDetections[clipNoiseIndexMask]
    noiseDetectedIndices = noiseDetectedIndices[clipNoiseIndexMask]
    peakDetections = peakDetections[clipPeakIndexMask]
    peakDetectedIndices = peakDetectedIndices[clipPeakIndexMask]
    print len(peakDetectedIndices)

    noiseDetectedIndicesShifted = noiseDetectedIndices+startPeakIndex-startNoiseIndex

    nRadioIndexBins = np.max(peakDetectedIndices)-np.min(peakDetectedIndices)

    histPeakIndices,peakIndexBins = np.histogram(peakDetectedIndices,bins=nRadioIndexBins)
    print peakIndexBins
    histNoiseIndices,noiseIndexBins = np.histogram(noiseDetectedIndices,bins=nRadioIndexBins)
    peakPhaseBins = indexToPhase(peakIndexBins[0:-1])


    radioIndexBins = np.arange(130.5,200,1)
    peakDist,distBinEdges = np.histogram(peakDetectedIndices,bins=radioIndexBins)
    print distBinEdges,peakDist
    noiseDist,distBinEdges = np.histogram(noiseDetectedIndicesShifted,bins=radioIndexBins)
    radioPhaseBins = indexToPhase(radioIndexBins)

    #plot the histograms in linear scale
    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.step(radioIndexBins[0:-1],noiseDist,'g',label='noise detections')
    ax.plot(radioPhaseBins,np.append(noiseDist,noiseDist[-1]),'g',drawstyle='steps-post',label='noise detections')
    ax.plot([radioPhaseBins[0],radioPhaseBins[-1]],np.repeat(np.mean(noiseDist),2),'g',drawstyle='steps-post',label='noise detections')
    #ax.step(radioIndexBins[0:-1],peakDist-noiseDist,'b',label='\'GRP detections\'')
    ax.plot(radioPhaseBins,np.append(peakDist,peakDist[-1]),'k',drawstyle='steps-post',label='GRP+noise detections')
    #ax.step(radioIndexBins[0:-1],peakDist,'k',label='GRP+noise detections')
    ax.legend(loc='best')
    ax.set_xlabel('radio phase')
    ax.set_title('Detection strength histograms')


    fig = plt.figure()
    ax = fig.add_subplot(111)
    #ax.step(radioIndexBins[0:-1],noiseDist,'g',label='noise detections')
    ax.plot(radioIndexBins,np.append(noiseDist,noiseDist[-1]),'g',drawstyle='steps-post',label='noise detections')
    ax.plot([radioIndexBins[0],radioIndexBins[-1]],np.repeat(np.mean(noiseDist),2),'g',drawstyle='steps-post',label='noise detections')
    #ax.step(radioIndexBins[0:-1],peakDist-noiseDist,'b',label='\'GRP detections\'')
    ax.plot(radioIndexBins,np.append(peakDist,peakDist[-1]),'k',drawstyle='steps-post',label='GRP+noise detections')
    #ax.step(radioIndexBins[0:-1],peakDist,'k',label='GRP+noise detections')
    ax.legend(loc='best')
    ax.set_xlabel('radio phase')
    ax.set_title('Detection strength histograms')
    #the probability that a radio detection at a certain radio strength is spurious is the ratio of total false detections in the test range to the total detections (true+false) in the P2 range
    spuriousProb = np.mean(noiseDist)/(1.*peakDist)
    spuriousProb[spuriousProb>1.]=1.
    grpProb = 1-spuriousProb

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(radioPhaseBins[0:-1],grpProb,'b',label='P(detection is GRP)')
    ax.plot(radioPhaseBins[0:-1],spuriousProb,'g',label='P(detection is noise)')
    ax.set_xlabel('radio arrival time')
    ax.legend(loc='upper right')

    #only include GRP that fall in the desired index range and strength range (threshold)
    radioPeakMask = np.array(clipPeakIndexMask)

    grpProbF = scipy.interpolate.interp1d(radioPhaseBins[0:-1],grpProb,kind='cubic')
    print 'interpolation range',radioPhaseBins[0],radioPhaseBins[-2]
    #grpProbAtIndices = grpProb[indices]
    nPoints = 100
    filledBins = np.linspace(radioPhaseBins[0],radioPhaseBins[-2],nPoints)
    fineSampledProb = grpProbF(filledBins)
    ax.plot(filledBins,fineSampledProb,'k--')


    return {'radioPeakMask':radioPeakMask,'radioPhaseBins':radioPhaseBins,'radioIndexBins':radioIndexBins,'grpProbFunc':grpProbF,'grpProbFuncRange':[radioPhaseBins[0],radioPhaseBins[-2]],'noiseDist':noiseDist,'peakDist':peakDist}




if __name__=='__main__':
    binsString = """0.130001     0.13137848   0.13254296   0.13372344   0.13481184
       0.1358914    0.13707388   0.13831236   0.13966668   0.14120796
          0.1428228    0.14459556   0.14662752   0.14906848   0.15211916
             0.1555064    0.160049     0.16608192   0.17604684   0.19213508
                0.222151     0.26999828   0.3452324    0.477951     0.75829176
                  23.311523 """
    radioStrengthBins = np.array([np.double(s) for s in binsString.split()])
    probsOfGRP(nRadioBins=15,strengthMin=.175)
    plt.show()
