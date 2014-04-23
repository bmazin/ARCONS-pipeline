import os, sys
from util import FileName
from util.ObsFile import ObsFile
from util.readDict import readDict
param = readDict()

INTERM_DIR = '.'
os.environ['INTERM_DIR'] = INTERM_DIR

run = 'PAL2012'
sunsetDate = '20121211'
obsSequence = '20121212-074700'
fn = FileName.FileName(run=run,date=sunsetDate,tstamp=obsSequence)
obsFile = ObsFile(fn.obs())
if len(sys.argv) > 1:
    obsFile.loadStandardCosmicMask()
    print " length of cosmicMask=",len(obsFile.cosmicMask)
firstSec = 0
integrationTime = 3
nPhoton = 0
for iRow in xrange(obsFile.nRow):
    for iCol in xrange(obsFile.nCol):
        gtpl = obsFile.getTimedPacketList(iRow,iCol,firstSec,integrationTime)
        nPhoton += len(gtpl['timestamps'])
print "nPhoton=",nPhoton
del obsFile
