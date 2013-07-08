import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
from util import hgPlot
from cosmic.Cosmic import Cosmic
import tables
import matplotlib.pyplot as plt
from hotpix import hotPixels
from scipy.optimize import curve_fit
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
            #print "seq=%s time=%12d contents=%3d"%(seq,ct,bc)
    
            run = 'PAL2012'
            sundownDate = '20121211'
            obsDate = '20121212'
            average = ct.mean()
            t0 = average-25
            t1 = t0+50
            plotfn = "cp-%05d-%5s-%5s-%5s-%5s-%010d"%(bc,run,sundownDate,obsDate,seq,t0)
            fn = FileName(run, sundownDate, obsDate+"-"+seq)
            cosmic = Cosmic(fn)
            dictionary = cosmic.fitExpon(t0, t1)
            plt.clf()
            hist = dictionary['timeHgValues']
            bins = np.arange(len(hist)) 
            plt.plot(bins, hist, label="event histogram")
            #mean = hist.mean()
            #plots the mean of the histogram
            #plt.vlines(mean, 0, 50, label='mean')
            

            xFit = np.linspace(0,len(hist),10*len(hist))
            pFit = dictionary['pFit']


            def funcExpon(x, a, b, c, d):
                retval = a*np.exp(-b*(x-d)) + c
                retval[x < d] = 0
                return retval
            def funcGauss(x, a, b, c):
                return a*np.exp(-(x-b)**2/(2.*c**2))

            pGaussFit = dictionary['pGaussFit']
            #pGaussFit, pcov = curve_fit(funcGauss, center, hist, 
             #                           p0=pGaussGuess)
       
            center = (bins[:-1]+bins[1:])/2
            histGaussFit = funcGauss(center,*pGaussFit)
            yFit = funcExpon(xFit,*pFit)
            plt.plot(xFit,yFit, label="exponential fit")
            plt.plot(center,histGaussFit, label="gaussian fit")
            plt.legend()
            plt.title(plotfn)
            plt.savefig(plotfn+".png")
