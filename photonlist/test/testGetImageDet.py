import time
from util.FileName import FileName
from util.utils import plotArray
from photonlist import photlist


def testGetImageDet(fileName=FileName(run='PAL2012',date='20121208',tstamp='20121209-120530').photonList(),
                    firstSec=0,integrationTime=-1,newMethod=True):
    
    plFile = photlist.PhotList(fileName)
    
    try:
        tic=time.time()
        image = plFile.getImageDet(firstSec=firstSec,integrationTime=integrationTime,newMethod=newMethod,
                                   wvlMin=4000,wvlMax=6000)
        tElapsed = time.time()-tic
    finally:
        #plFile.close()
        print 'Deleting file instance'
        del plFile
    plotArray(image)
    print 'Done, time taken (s): ',tElapsed
    return image