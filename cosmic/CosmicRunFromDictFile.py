import sys, logging, pickle
from util.ObsFile import ObsFile
from cosmic.Cosmic import Cosmic
from util.readDict import readDict
from util import FileName
import matplotlib.pyplot as plt
import numpy as np
class CosmicRunFromDictFile:
    def __init__(self, dictFile, cosmicDictFile):
        # define which files to process
        self.dictFile = dictFile
        self.param = readDict()
        self.param.readFromFile(dictFile)

        # define cosmic ray finding parameters
        self.cosmicDictFile = cosmicDictFile
        self.s = readDict()
        self.s.readFromFile(cosmicDictFile)

    def run(self):
        for obs in self.param['obsSequence']:
            fn = FileName.FileName(run=self.param['run'],
                                   date=self.param['sunsetDate'],
                                   tstamp = obs)


            cosmic = Cosmic(fn,
                            loggingLevel=logging.INFO,
                            beginTime=self.s['beginTime'],
                            endTime=self.s['endTime'])


            fc = cosmic.findCosmics(stride=int(self.s['stride']), 
                                    threshold=int(self.s['threshold']), 
                                    populationMax=int(self.s['populationMax']),
                                    nSigma=float(self.s['nSigma']),
                                    writeCosmicMask=True)

            cosmic = Cosmic(fn,
                            loggingLevel=logging.INFO,
                            beginTime=self.s['beginTime'],
                            endTime=self.s['endTime'],
                            applyCosmicMask=True)

            fcm = cosmic.findCosmics(stride=int(self.s['stride']), 
                                     threshold=int(self.s['threshold']), 
                                     populationMax=int(self.s['populationMax']),
                                     nSigma=float(self.s['nSigma']),
                                     writeCosmicMask=False)
            
            p = {
                'populationHg':fc['populationHg'][0],
                'populationHgM':fcm['populationHg'][0],
                'pps':fc['pps'],
                'ppsM':fcm['pps'],
                'ppmMasked':fc['ppmMasked'],
                'ppsTime':fc['ppsTime'],
                'beginTime':cosmic.beginTime,
                'endTime':cosmic.endTime
                }

            outfileName = "cosmic-summary-"+obs+".pkl"
            outfile = open(outfileName, "wb")
            pickle.dump(p,outfile)
            outfile.close()

            
    def summarize(self):
        print "summarize"
        nSum = 0
        for obs in self.param['obsSequence']:
            pFileName = "cosmic-summary-"+obs+".pkl"
            pFile = open(pFileName, "rb")
            p = pickle.load(pFile)
            pFile.close()
            stride = int(1/p['ppsTime'])

            if nSum > 0:
                popAllSum += p['populationHg']
                popMaskedSum += p['populationHgM']
                ppmMaskedSum += p['ppmMasked']
                pps = p['pps']
                ppsSum += np.sum(pps.reshape(-1,stride),axis=0)
                ppsMasked = p['ppsM']
                ppsMaskedSum += np.sum(ppsMasked.reshape(-1,stride),axis=0)
            else:
                popAllSum = p['populationHg']
                popMaskedSum = p['populationHgM']
                ppmMaskedSum = p['ppmMasked']
                pps = p['pps']
                ppsSum = np.sum(pps.reshape(-1,stride),axis=0)
                ppsMasked = p['ppsM']
                ppsMaskedSum = np.sum(ppsMasked.reshape(-1,stride),axis=0)
            nSum += 1

        ppmMasked = ppmMaskedSum/float(nSum)
        plt.figure()
        n = len(ppsSum)
        print "in plot pps.png n = ",n
        xTime = np.arange(n)/float(n)
        plt.plot(xTime, ppsSum,label="all data")
        plt.plot(xTime, ppsMaskedSum, label="cosmic masked")
        plt.legend(loc="best").get_frame().set_alpha(0.5)
        plt.title("%s ppmMasked=%d"%(self.dictFile,ppmMasked))
        plt.xlabel("time inside each second")
        plt.ylabel("count rate")
        plt.savefig("pps.png")
    
        plt.figure()
        plt.plot(popAllSum, "b-", drawstyle="steps-mid", label="all")
        plt.plot(popMaskedSum, "r+", drawstyle="steps-mid", label="cosmics masked")

        plt.xscale("symlog", linthreshx=0.9)
        plt.yscale("symlog", linthreshy=0.1)
        plt.xlim(-0.1, self.s['populationMax'])
        plt.ylim(ymin=-0.1)
        plt.legend(loc="center right").get_frame().set_alpha(0.5)
        plt.title("%s ppmMasked=%d"%(self.dictFile,ppmMasked))
        plt.xlabel("Population -- number of photons in a time tick")
        plt.ylabel("dN/dPopulation")
        plt.savefig("populationHg.png")

               

if __name__ == '__main__':
    try:
        dictFile = sys.argv[1]
        cosmicDictFile = sys.argv[2]
    except IndexError:
        print "Usage:  CosmicRunFromDiceFile.py dictFile cosmicDictFile"
        print " to run cosmic ray finding code"
        print " "
        print "or"
        print "Usage:  CosmicRunFromDiceFile.py dictFile cosmicDictFile 1"
        print " to summarize results of cosmic ray finding code"
        exit()
    crfdf = CosmicRunFromDictFile(dictFile, cosmicDictFile)
    if len(sys.argv) == 3:
        crfdf.run()
    else:
        crfdf.summarize()
