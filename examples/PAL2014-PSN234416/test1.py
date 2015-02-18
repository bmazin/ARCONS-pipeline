from util.ObsFile import ObsFile
from util.FileName import FileName
run = "PAL2014"
date = "20141022"
timeStamp = '20141023-033821'
fn = FileName(run,date,timeStamp)

of = ObsFile(fn.obs())
of.loadBestWvlCalFile()
print "wvlCalFileName=",of.wvlCalFileName
of.loadBeammapFile(fn.beammap())
fn2 = FileName(run,date,"")
of.loadFlatCalFile(fn2.flatSoln())
row = 4
col = 4
firstSec = 72
integrationTime = 73
spec = of.getPixelSpectrum(row,col,firstSec,integrationTime)
print "spec=",spec

