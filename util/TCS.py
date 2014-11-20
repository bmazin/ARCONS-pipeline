import numpy as np
from util import FileName
from util import ObsFile
class TCS:
    """
    Access the TCS logs to get the ra,dec motions for a range of observations
    """
    def __init__(self,run,date):
        """
        The run ("LICK2014") and date ("20141020") define the entire log file
        """
        tempFileName = FileName.FileName(run,date)
        tcsPath = tempFileName.tcslog()
        with open(tcsPath,'r') as f:
            tcsLines = f.readlines()
        self.time = np.zeros(len(tcsLines))
        self.dra = np.zeros(len(tcsLines))
        self.ddec = np.zeros(len(tcsLines))
        self.raOffset = np.zeros(len(tcsLines)+1)
        self.decOffset = np.zeros(len(tcsLines)+1)
        
        for i,line in enumerate(tcsLines):
            lineList = line.split()
            self.time[i] = float(lineList[0])
            self.dra[i] = float(lineList[3])
            self.ddec[i] = float(lineList[4])

    def select(self,obsBegin,obsEnd):
        """
        Select all the moves that occur from the beginning of obsBegin
        to the end of obsEnd, given by the ObsFile

        Return a dictionary of "time" (UNIX time in seconds) ,"dra", and "ddec"
        (in arcsec)
        """
        t0 = obsBegin.getFromHeader('unixtime')
        t1 = obsEnd.getFromHeader('unixtime')+obsEnd.getFromHeader('exptime')
        b = (t0 < self.time) & (self.time < t1)
        times = self.time[b]
        dras = self.dra[b]
        ddecs = self.ddec[b]
        raOffsets = np.zeros(len(times)+1)
        decOffsets = np.zeros(len(times)+1)
        for i in range(len(times)):
            raOffsets[i+1] = raOffsets[i]+dras[i]
            decOffsets[i+1] = decOffsets[i]+ddecs[i]
        return {"time":times, "dra":dras, "ddec":ddecs,
                "raOffset":raOffsets,"decOffset":decOffsets}

# Simple test
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    run = 'PAL2014'
    date = '20141020'
    tcs = TCS(run,date)
    dt0 = '20141021-045212'
    dt1 = '20141021-051726'
    of0 = ObsFile.ObsFile(FileName.FileName(run,date,dt0).obs())
    of1 = ObsFile.ObsFile(FileName.FileName(run,date,dt1).obs())    
    tcsDict = tcs.select(of0,of1)
    del of0
    del of1
    dx = tcsDict['dra']
    dy = tcsDict['ddec']
    x = 0
    y = 0
    xl = [x]
    yl =[y]
    for i in range(len(dx)):
        x += dx[i]
        y += dy[i]
        xl.append(x)
        yl.append(y)
        plt.text(x,y,str(i+1),
                 horizontalalignment="center",
                 verticalalignment="center")
    plt.plot(xl,yl)
    plt.xlabel("delta ra (arcsec)")
    plt.ylabel("delta dec (arcsec)")
    plt.xlim(min(xl)-1,max(xl)+1)
    plt.ylim(min(yl)-1,max(yl)+1)
    plt.show()
