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

            if True:
                cosmic = Cosmic(fn,
                                loggingLevel=logging.INFO,
                                beginTime=self.s['beginTime'],
                                endTime=self.s['endTime'])


                fc = cosmic.findCosmics(stride=int(self.s['stride']), 
                                        threshold=int(self.s['threshold']), 
                                        populationMax=int(self.s['populationMax']),
                                        nSigma=float(self.s['nSigma']),
                                        writeCosmicMask=True)

                outfileName = "cosmic-all"+obs+".pkl"
                outfile = open(outfileName, "wb")
                pickle.dump(fc['populationHg'],outfile)
                pickle.dump(cosmic.beginTime,outfile)
                pickle.dump(cosmic.endTime,outfile)
                pickle.dump(fc['ppmMasked'],outfile)
                pickle.dump(fc['ppsTime'],outfile)
                pickle.dump(fc['pps'],outfile)
                outfile.close()

            cosmic = Cosmic(fn,
                            loggingLevel=logging.INFO,
                            beginTime=self.s['beginTime'],
                            endTime=self.s['endTime'],
                            applyCosmicMask=True)

            fc = cosmic.findCosmics(stride=int(self.s['stride']), 
                                    threshold=int(self.s['threshold']), 
                                    populationMax=int(self.s['populationMax']),
                                    nSigma=float(self.s['nSigma']),
                                    writeCosmicMask=False)
            
            outfileName = "cosmic-masked"+obs+".pkl"
            outfile = open(outfileName, "wb")
            pickle.dump(fc['populationHg'],outfile)
            pickle.dump(cosmic.beginTime,outfile)
            pickle.dump(cosmic.endTime,outfile)
            pickle.dump(fc['ppmMasked'],outfile)
            pickle.dump(fc['ppsTime'],outfile)
            pickle.dump(fc['pps'],outfile)
            outfile.close()

            
    def summarize(self):
        print "summarize"
        popAllSum  = None
        popMaskedSum = None
        ppsTimeMasked = {}
        ppsMasked = {}
        ppsTimeAll = {}
        ppsAll = {}
        for obs in self.param['obsSequence']:
            fn = FileName.FileName(run=self.param['run'],
                                   date=self.param['sunsetDate'],
                                   tstamp = obs)

            pFileName = "cosmic-all"+obs+".pkl"
            pFile = open(pFileName, "rb")
            popAll = pickle.load(pFile)
            beginTimeAll = pickle.load(pFile)
            endTimeAll = pickle.load(pFile)
            ppmMaskedAll = pickle.load(pFile)
            ppsTimeAll[obs] = pickle.load(pFile)
            ppsAll[obs] = pickle.load(pFile)

            pFileName = "cosmic-masked"+obs+".pkl"
            pFile = open(pFileName, "rb")
            popMasked = pickle.load(pFile)
            beginTimeMasked = pickle.load(pFile)
            endTimeMasked = pickle.load(pFile)
            ppmMaskedMasked = pickle.load(pFile)
            ppsTimeMasked[obs] = pickle.load(pFile)
            ppsMasked[obs] = pickle.load(pFile)

            cfn = fn.cosmicMask()
            table = ObsFile.readCosmicIntervalFromFile(cfn)

            if popAllSum:
                popAllSum += popAll[0]
                popMaskedSum += popMasked[0]
                ppmMaskedSum += ppmMaskedMasked
            else:
                popAllSum = popAll[0]
                popMaskedSum = popMasked[0]
                ppmMaskedSum = ppmMaskedMasked

        plt.figure()
        plt.plot(popAllSum, "b-", drawstyle="steps-mid", label="all")
        plt.plot(popMaskedSum, "r+", drawstyle="steps-mid", label="cosmics masked")

        plt.xscale("symlog", linthreshx=0.9)
        plt.yscale("symlog", linthreshy=0.1)
        plt.xlim(-0.1, self.s['populationMax'])
        plt.ylim(ymin=-0.1)
        plt.legend(loc="center right").get_frame().set_alpha(0.5)
        plt.title("%s ppmMasked=%d"%(self.dictFile,ppmMaskedSum))
        plt.xlabel("Population -- number of photons in a time bin")
        plt.ylabel("dN/dPopulation")

        plt.savefig("summary.png")

        plt.figure()
        for obs in self.param['obsSequence']:
            time = np.arange(len(ppsAll[obs]))*ppsTimeAll[obs]
            plt.plot(time,ppsAll[obs]-ppsMasked[obs])
        plt.savefig("pps.png")    
               

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
