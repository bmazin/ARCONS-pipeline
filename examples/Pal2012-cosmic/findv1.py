
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
threshold = 40
populationMax=2000
ySum = np.zeros(populationMax)
frameSum = 'none'

for seq in seq5:
    print "seq=",seq
    fn = FileName(run, sundownDate, obsDate+"-"+seq)
    cosmic = Cosmic(fn, endTime='exptime')
    fc = cosmic.findCosmics(stride, threshold, populationMax)
    if frameSum == 'none':
        frameSum = fc['frameSum']
    else:
        frameSum += fc['frameSum']
    outfile = open("cosmicTimeList-"+seq+".pkl", "wb")
    pickle.dump(fc['cosmicTimeList'],outfile)
    pickle.dump(fc['binContents'],outfile)
    outfile.close()
    populationHg = fc['populationHg']
    yPlot = populationHg[0].astype(np.float)
    ySum += yPlot
    norm = yPlot.sum()
    yPlot /= norm
    yPlot[yPlot==0] = 1e-3/norm
    xPlot = 0.5*(populationHg[1][1:]+populationHg[1][:-1])
    plt.clf()
    plt.plot(xPlot,yPlot,'-')
    plt.yscale('log')
    plt.ylim(ymin=0.5/norm)
    plt.xlim(1,populationMax)
    plt.xscale('log')
    #poisson = []
    #xValues = populationHg[1][1:]-0.5
    #for x in xValues:
    #    prob = (mean**x)*math.exp(-mean)/math.factorial(x)
    #    poisson.append(prob)
    #    pError = np.sqrt(poisson/sum(populationHg[0]))
    #    plt.errorbar(xValues, poisson, pError, label=
    plt.title("%s %s %s %s s=%d"%(run,sundownDate,obsDate,seq,stride))
    plt.savefig(seq+".png")
    del cosmic
    #if seq == '113212':
    #    print "break now"
    #    break

plt.clf()
ySum[ySum==0] = 1e-3
plt.plot(xPlot,ySum,'-')
plt.yscale('log')
plt.xlim(1,populationMax)
plt.xscale('log')
plt.ylim(ymin=0.5)
plt.title("%s %s %s-%s stride=%d"%(run,sundownDate,seq5[0],seq5[-1],stride))
plt.xlabel("number of simultaneous photons in array")
plt.ylabel("dN/dx")
plotName = "sum-%s-%s-%s-stride-%d"%(run,seq5[0],seq5[-1],stride)
plt.savefig(plotName)

frameFile = open("frameSum.pkl","wb")
pickle.dump(frameSum,frameFile)
frameFile.close()
