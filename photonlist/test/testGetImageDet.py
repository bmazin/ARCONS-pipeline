import time
from util.utils import plotArray
from photonlist import photlist


def testGetImageDet(fileName=None,firstSec=0,integrationTime=-1,newMethod=True):
    
    plFile = photlist.PhotList(fileName)
    
    try:
        tic=time.time()
        image = plFile.getImageDet(firstSec=firstSec,integrationTime=integrationTime,newMethod=newMethod)
        tElapsed = time.time()-tic
    finally:
        #plFile.close()
        print 'Deleting file instance'
        del plFile
    plotArray(image)
    print 'Done, time taken (s): ',tElapsed
    return image