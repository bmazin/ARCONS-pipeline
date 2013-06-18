import numpy as np

#This is used to combine multiple Image Stack npz files into one large file.

npzfile = np.load('/Scratch/dataProcessing/SDSS_J0926/AllData/Dec8ShortIntImageStackRed.npz')
npzfile2 = np.load('/Scratch/dataProcessing/SDSS_J0926/AllData/Dec10SIImageStackRednewCal.npz')
npzfile3 = np.load('/Scratch/dataProcessing/SDSS_J0926/AllData/Dec11SIImageStackRed.npz')
#npzfile.files

As=np.array(npzfile['stack'])
Aj=np.array(npzfile['jd'])
Bs=np.array(npzfile2['stack'])
Bj=np.array(npzfile2['jd'])
Cs=np.array(npzfile2['stack'])
Cj=np.array(npzfile2['jd'])

S = np.concatenate((As,Bs),axis = 2)
J = np.concatenate((Aj,Bj))

np.savez('/Scratch/dataProcessing/SDSS_J0926/AllData/AllDataSIImageStackRednewCal.npz',stack=S,jd=J)
