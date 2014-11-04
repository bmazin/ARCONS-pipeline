import numpy as np


labels='iBTOA   pulseNumber BTOA    Noise_Offset    Noise_RMS   Max  Mean Index    TestMax  TestMean TestIndex'
labels = np.array(labels.split())
path = '/Scratch/dataProcessing/crabData2/'

p2=np.loadtxt(path+'giantPulseList_P2_3sigma_indices.txt',skiprows=1)
p1=np.loadtxt(path+'giantPulseList_P1_3sigma_indices.txt',skiprows=1)

p1 = p1[p1[:,np.argmax(labels=='Max')]>=.155]
p2 = p2[p2[:,np.argmax(labels=='Max')]>0.]

print len(p1),len(p2)

giantDictP1 = dict()
giantDictP2 = dict()
for iLabel,label in enumerate(labels):
    giantDictP1[label] = p1[:,np.argmax(label==labels)]
    giantDictP2[label] = p2[:,np.argmax(label==labels)]


pulseNumbersP1 = giantDictP1['pulseNumber']
pulseNumbersP2 = giantDictP2['pulseNumber']
overlap = []
for pn in pulseNumbersP1:
    if pn in pulseNumbersP2:
        overlap.append(pn)
overlap  = np.array(overlap)
print len(pulseNumbersP1),len(pulseNumbersP2),len(overlap)
np.savez('overlapP1.npz',overlap=overlap,p1=p1,p2=p2)

