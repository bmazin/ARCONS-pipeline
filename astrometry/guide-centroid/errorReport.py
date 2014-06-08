import pyfits 
from functions import *
import os

fileList = []
for nfile in os.listdir('./cal/'):
    try:
        if str(nfile[6:]) == '_offCal_rotCal.fits':
            fileList.append(nfile)
    except:
        pass

error = []
tol = 2000
count = 0
for nfile in fileList:
    imageList = pyfits.open('./cal/'+nfile)
    header = imageList[0].header
    if header['CALERR'] > 2000:
        error.append(fileList[count][0:6])
    imageList.close()
    count += 1
print 'total files = %s and %s of them are bad!' %(len(fileList),len(error))
print error