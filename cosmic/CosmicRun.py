import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
from util import hgPlot
from cosmic.Cosmic import Cosmic
import tables
from scipy.optimize import curve_fit
from hotpix import hotPixels
import pickle
from interval import interval, inf, imath
from cosmic import tsBinner
import os
import sys
import logging
class CosmicRun:
    
    def __init__(self, path):
        print "begin path=",path
        os.chdir(path)
        file = open("settings.par", 'r')
        self.s = {}
        for line in file:
            temp = line.split("=")
            if len(temp) > 1:
                self.s[temp[0].strip()] = temp[1].strip()
                print temp[0].strip(), '=', temp[1].strip()
        file.close()

    def findv1(self):
        populationMax=2000
        ySum = np.zeros(populationMax)
        frameSum = 'none'
        seq5 = self.s['seq5'].split()
        for seq in seq5:
            print "seq=",seq
            outfileName = "cosmicTimeList-"+seq+".pkl"
            if not os.path.exists(outfileName):
                fn = FileName(self.s['run'], 
                              self.s['sundownDate'], 
                              self.s['obsDate']+"-"+str(seq))
                cosmic = Cosmic(fn, 
                                beginTime=self.s['beginTime'],
                                endTime=self.s['endTime'],
                                loggingLevel = logging.INFO)

                fc = cosmic.findCosmics(stride=int(self.s['stride']), 
                                        threshold=int(self.s['threshold']), 
                                        populationMax=populationMax,
                                        nSigma=float(self.s['nSigma']))
                outfile = open(outfileName, "wb")
                pickle.dump(fc['cosmicTimeList'],outfile)
                pickle.dump(fc['binContents'],outfile)
                outfile.close()
                cfn = "cosmicMask-%s.h5"%seq
                ObsFile.writeCosmicIntervalToFile(fc['interval'],1.0e6, cfn,
                                                  self.s['beginTime'],
                                                  self.s['endTime'],
                                                  int(self.s['stride']),
                                                  int(self.s['threshold']),
                                                  float(self.s['nSigma']),
                                                  populationMax)

                del cosmic

    def makemovie1(self):
        run = self.s['run']
        sundownDate = self.s['sundownDate']
        obsDate = self.s['obsDate']
        stride = int(self.s['stride'])
        seq5 = self.s['seq5'].split()
        for seq in seq5:
            inFile = open("cosmicTimeList-%s.pkl"%(seq),"rb")
            cosmicTimeList = pickle.load(inFile)
            binContents = pickle.load(inFile)
            cfn = "cosmicMask-%s.h5"%seq
            intervals = ObsFile.readCosmicIntervalFromFile(cfn)
            for interval in intervals:
                print "interval=",interval
                fn = FileName(run, sundownDate,obsDate+"-"+seq)
                obsFile = ObsFile(fn.obs())
                obsFile.loadTimeAdjustmentFile(fn.timeAdjustments())
                i0=interval[0]
                i1=interval[1]
                intervalTime = i1-i0
                dt = intervalTime/2
                beginTime = max(0,i0-0.000200)
                endTime = beginTime + 0.001
                integrationTime = endTime-beginTime
                nBins = int(np.round(obsFile.ticksPerSec*(endTime-beginTime)+1))
                timeHgValues = np.zeros(nBins, dtype=np.int64)
                ymax = sys.float_info.max/100.0
                for iRow in range(obsFile.nRow):
                    for iCol in range(obsFile.nCol):
                        gtpl = obsFile.getTimedPacketList(iRow,iCol,
                                                          beginTime,integrationTime)
                        ts = (gtpl['timestamps'] - beginTime)*obsFile.ticksPerSec
                        ts64 = np.round(ts).astype(np.uint64)
                        tsBinner.tsBinner(ts64, timeHgValues)
                plt.clf()
                plt.plot(timeHgValues, label="data")
                x0 = (i0-beginTime)*obsFile.ticksPerSec
                x1 = (i1-beginTime)*obsFile.ticksPerSec
                plt.fill_between((x0,x1),(0,0), (ymax,ymax), alpha=0.2, color='red')   
                plt.yscale("symlog",linthreshy=0.9)
                plt.xlim(0,1000)
                plt.ylim(-0.1,300)
                tick0 = int(np.round(i0*obsFile.ticksPerSec))
                plotfn = "cp-%05d-%s-%s-%s-%09d"%(timeHgValues.sum(),run,obsDate,seq,tick0)
                plt.title(plotfn)
                plt.legend()
                plt.savefig(plotfn+".png")
                plt.xlabel("nSigma=%d stride=%d threshold=%d"%(int(self.s['nSigma']),int(self.s['stride']),int(self.s['threshold'])))
                print "plotfn=",plotfn
                

        os.system("convert -delay 0 `ls -r cp*png` cp.gif")




if __name__ == '__main__':
    if len(sys.argv) >1:
        path = sys.argv[1]
    else:
        path = "."
    cosmicRun = CosmicRun(path)
    cosmicRun.findv1()
    print "now call makemovie1"
    cosmicRun.makemovie1()
    print "glorious success"
