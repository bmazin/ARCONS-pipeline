import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
from util import hgPlot
from cosmic.Cosmic import Cosmic
import tables
import matplotlib.pyplot as plt
from hotpix import hotPixels
import pickle

run = 'PAL2012'
sundownDate = '20121211'
obsDate = '20121212'
# December 11
# Final sequence, toward end of run, thin high clouds at around 12:50, moved to confirm position at '122234', also at '130752' at 125s. (16,15)
seq5 = ['112709', '113212', '113714', '114216', '114718', '115220', '115722', '120224', '120727', '121229', '121732', '122234', '122736', '123238', '123740', '124242', '124744', '125246', '125748', '130250', '130752', '131254', '131756', '132258', '132800', '133303']

stride = 10
threshold = 600
nAboveThreshold = 0

for seq in seq5:
    inFile = open("cosmicTimeList-%s.pkl"%(seq),"rb")
    cosmicTimeList = pickle.load(inFile)
    binContents = pickle.load(inFile)
    for i in range(cosmicTimeList.size):
        bc = binContents[i]
        if bc > threshold:
            nAboveThreshold += 1
            ct = cosmicTimeList[i]
            print "seq=%s time=%12d contents=%3d"%(seq,ct,bc)
print "nAboveThreshold=",nAboveThreshold
