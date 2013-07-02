import sys
from util.FileName import FileName
from cosmic.Cosmic import Cosmic
from util import utils
import pickle
import matplotlib as mpl
import numpy as np
run = 'PAL2012'
sundownDate = '20121211'
obsDate = '20121212'

args = sys.argv
seq = args[1]
tick0 = np.uint64(args[2])
sec0 = tick0/1000000
print "seq=%s tick0=%d sec0=%d"%(seq,tick0,sec0)

fn = FileName(run, sundownDate, obsDate+"-"+seq)
cosmic = Cosmic(fn, endTime='exptime')

frameFile = open("frameSum.pkl","rb")
frameSum = pickle.load(frameFile)
frameFile.close()
print "frameSum.shape=",frameSum.shape

frames = []
titles = []
t0 = tick0-np.uint64(10)
t1 = tick0+np.uint64(60)

def makeTitle(tick,t0,t1):
    l = ['.']*(t1-t0)
    l[np.uint64(tick)-t0] = "*"
    title = "tick=%10d usec %s"%(tick,"".join(l))
    return title
listOfPixelsToMark = []
for tick in range (t1-t0):
    listOfPixelsToMark.append([])

for iRow in range(cosmic.file.nRow):
    for iCol in range(cosmic.file.nCol):
        gtpl = cosmic.file.getTimedPacketList(iRow,iCol,sec0,1)
        timestamps = gtpl['timestamps']
        timestamps *= cosmic.file.ticksPerSec
        ts64 = timestamps.astype(np.uint64)
        for ts in ts64:
            tindex = ts-t0
            try:
                listOfPixelsToMark[tindex].append((iRow,iCol))
            except IndexError:
                pass
for tick in range(t0,t1):
    frames.append(frameSum)
    title = makeTitle(tick,t0,t1)
    titles.append(title)

mfn0 = "m-%s-%s-%s-%s-%010d-%010d-i.gif"%(run,sundownDate,obsDate,seq,t0,t1)
utils.makeMovie(frames, titles, outName=mfn0, delay=0.1, colormap=mpl.cm.gray,
                listOfPixelsToMark=listOfPixelsToMark,
                pixelMarkColor='red')

for i in range(len(listOfPixelsToMark)-1):
    listOfPixelsToMark[i+1].extend(listOfPixelsToMark[i])

mfn1 = "m-%s-%s-%s-%s-%010d-%010d-a.gif"%(run,sundownDate,obsDate,seq,t0,t1)
utils.makeMovie(frames, titles, outName=mfn1, delay=0.1, colormap=mpl.cm.gray,
                listOfPixelsToMark=listOfPixelsToMark,
                pixelMarkColor='green')

if __name__ == "__main__":
    import sys
    
