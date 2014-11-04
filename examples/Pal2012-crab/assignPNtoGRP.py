import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.stdout.flush()
np.set_printoptions(precision=11)
path = '/Scratch/dataProcessing/crabData/'
path2 = '/Scratch/dataProcessing/crabData2/'

pulseLabel = 1
labels='BTOA Noise_Offset Noise_RMS Max Mean Index TestMax TestMean TestIndex'
labels = np.array(labels.split())
radioGiantData = np.loadtxt(path2+'radio/Giant_List_P{}_BTOA_Flux_Sorted'.format(pulseLabel),usecols=range(1,len(labels)+1))
giantDict = dict()
print np.shape(radioGiantData)
for iLabel,label in enumerate(labels):
    labelIdx = np.argmax(label==labels)
    giantDict[label] = radioGiantData[:,labelIdx]
radioGiantBTOA = giantDict['BTOA']#radioGiantData[:,0]
opticalData = np.load(path+'btoa20121211wave.npz')#has timeAdjustFile incorporated


opticalBTOA = opticalData['btoaList']
sortIdx = np.argsort(opticalBTOA)
opticalBTOA = opticalBTOA[sortIdx]
opticalCumulativePhase = opticalData['btoaPhase'][sortIdx]

opticalPulseNumbers = np.array(opticalCumulativePhase,dtype=np.int)
    
pulsarPeriod = 33e-3/(3600*24) #days, approximately
    #pulseStartIdx = np.searchsorted(opticalPulseNumbers,pulseIdx+idxOffset)

print 'iRadioBTOA radioBTOA opticalPulseNumber radioMax radioMean'
radioPulseNumbers = np.zeros(len(radioGiantBTOA),dtype=np.uint64)
periodThreshold=.5
outFile = open(os.path.join(path2,'giantPulseList_P{}_3sigma_indices.txt'.format(pulseLabel)),mode='w')
outFile.write('iBTOA\tpulseNumber')
for label in labels:#everything after BTOA
    outFile.write('\t'+label)
outFile.write('\n')
for iBTOA,radioBTOA in enumerate(radioGiantBTOA):
    idx = np.searchsorted(opticalBTOA,radioBTOA)
    pulseNumber = opticalPulseNumbers[idx]
    #check that the nearest optical photon and the radioBTOA
    # aren't more than half a pulse period apart
    # indicating that we have less than half (probably none) of the optical GRP
    periodsOff = (opticalBTOA[idx]-radioBTOA)/pulsarPeriod
    if (np.abs(periodsOff) < periodThreshold):
        print 'matched pulse',pulseNumber,'by',periodsOff
        radioPulseNumbers[iBTOA]=pulseNumber
        outFile.write(str(iBTOA)+'\t'+'\t'+str(pulseNumber))
        remainingRow = ''.join(['\t'+str(cell) for cell in radioGiantData[iBTOA,:]])
        outFile.write(remainingRow)
        outFile.write('\n')
    else:
        print 'missed pulse',radioBTOA

outFile.close()
print 'done'

