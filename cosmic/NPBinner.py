import numpy as np
class NPBinner:
    @staticmethod
    def binnerSlow(timeStamps,bins):
        for timeStamp in timeStamps:
            bins[timeStamp] += 1
        return bins

    @staticmethod
    def binnerNotAsSlow(timeStamps,nBins):
        hg = np.histogram(timeStamps,np.arange(nBins))
        return hg[0]

    @staticmethod
    def binner(timeStamps, nBins):
        contents = np.bincount(timeStamps, minlength=nBins)
        return contents
