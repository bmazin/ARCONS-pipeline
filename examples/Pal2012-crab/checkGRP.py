import numpy as np
import matplotlib.pyplot as plt
import sys
import os
sys.stdout.flush()
np.set_printoptions(precision=11)
path = '/Scratch/dataProcessing/crabData/'
path2 = '/Scratch/dataProcessing/crabData2/'

labels='BTOA Noise_Offset Noise_RMS Max_P2 Mean_P2 Index_P2 TestMax_P2 TestMean_P2 TestIndex_P2'
labels = np.array(labels.split())
radioGiantData = np.loadtxt(path2+'radio/Giant_List_P2_BTOA_Flux_Sorted',usecols=range(1,len(labels)+1))

radioGiantData2 = np.loadtxt(path+'radio/GiantPulses_BTOA_Flux_P2',usecols=range(1,len(labels)+1),skiprows=1)

print np.shape(radioGiantData),np.shape(radioGiantData2)
radioGiantData = radioGiantData[(radioGiantData[:,labels=='Max_P2']>=1.).reshape(-1)]
radioGiantData2 = radioGiantData2[(radioGiantData2[:,labels=='Max_P2']>=1.).reshape(-1)]
print np.shape(radioGiantData),np.shape(radioGiantData2)
giantDict = dict()
giantDict2 = dict()
for iLabel,label in enumerate(labels):
    labelIdx = np.argmax(label==labels)
    print iLabel,label,labelIdx
    giantDict[label] = radioGiantData[:,labelIdx]
    giantDict2[label] = radioGiantData2[:,labelIdx]

radioGiantBTOA = giantDict['BTOA']#radioGiantData[:,0]
radioGiantBTOA2 = giantDict2['BTOA']#radioGiantData[:,0]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(giantDict['BTOA'],giantDict['Max_P2'],'b.')
ax.plot(giantDict2['BTOA'],giantDict2['Max_P2'],'g.')
print 'done'
plt.show()

