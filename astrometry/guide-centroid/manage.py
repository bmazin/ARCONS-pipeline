from FitsAnalysis import convert,StarCalibration
from catalog import queryVizier,queryFitsImage
import os
import warnings
from radec import radec
from functions import *


#ignore the warning caused by astropy
warnings.filterwarnings("ignore")

#This specifies the center of fits images retrieved from the data base
#Though it is possible to specify by name, it is a good idea to use 'RA DEC' in degrees to avoid erros
#pos = '104.9566125,14.2341555'
#pos = 'corot18b'
#RA,DEC = convert([['6:59:49.587','+14:14:02.96']])[0]
#print RA,DEC
#pos = '%s,%s' %(RA,DEC)

#PSR 0656+14
pos = '104.95,14.24'

#source of data:'USNO-B1.0' or '2MASS' are usually enough. For full list: http://cdsarc.u-strasbg.fr/viz-bin/vizHelp?cats/U.htx
source = 'USNO-B1.0'
#source = '2MASS'

#name of the saved files
tfitsTable = 'test.fits'
tfitsImage = 'test_image.fits'

#if manCat=True, manual catalog will be used instead of vizier
#if semiManCat=True, stars will be added on top of the vizier catalog stars
#stars appended in both cases are specified in manCatFile
#notice that manCat and semiManCat can't both be true at the same time
manCat = False
semiManCat = True
manCatFile = 'manCat.cat'

calHeight = 3

#saving directory of all the calibrated files in relative path
caldir = './cal/'
#directory of fits images to be calibrated, put all the files here
fdir = './origin/'
sedir = './config/'
#the distoriton parameter file
paramFile = None

#if manual = False, the program will use sextractor to find source and match the correponding stars in the images
#also make sure the ./origin/ folder has appropriate sextractor parameters files and parameters
manual = False

#if calibrate is True, all the files that are calibrated will be used as data points to calculate distortion parameters
calibrate = False


#next, if automatic calibration is chosen, it is best to first manually correct the reference pixel coordinate on the header. This greatly increases the chances of calibrating.
refFix = True
CRVAL1 = 104.94975
CRVAL2 = 14.2423833333333


'''
-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------

Input Ends Here

-----------------------------------------------------------------------------------------------------------
-----------------------------------------------------------------------------------------------------------
'''

#it will overwrite any existing files with the same names, the file 'test_vo.xml' is not important and can be ignored
#queryVizier(tfitsTable,source=source,pos=pos)
#queryFitsImage(tfitsImage,'test_vo.xml',pos=pos)
if manCat and semiManCat:
    raise ValueError, 'Manual catalog and semi-manual catalog cannot be True all at once!'
elif manCat:
    catOption = 'full'
elif semiManCat:
    catOption = 'semi'
else:
    catOption = None
    

#perform linear and polynomial calibration to each file in dir specified
for fitsImage in os.listdir(fdir):
    
    #I am separation lines  
    print '--------------------------------------------------------------------------'
    print '--------------------------------------------------------------------------'
    print '> Calibrating %s...' %(fitsImage)
    
    #fix reference value if refFix is True
    if refFix:
        updateHeader(fdir+fitsImage,'CRVAL1',CRVAL1)
        updateHeader(fdir+fitsImage,'CRVAL2',CRVAL2) 
    try:
        cal = StarCalibration(fitsImage,tfitsTable,tfitsImage,manual,paramFile=paramFile,caldir=caldir,fdir=fdir,sedir=sedir,height=3,manCat=catOption,manCatFile=manCatFile)
        cal.linCal()
    
        if paramFile != None:
            distHeaderUpdate(caldir+fitsImage[:-5]+'_offCal_rotCal.fits',caldir+fitsImage[:-5]+'_allCal.fits',paramFile)
            #cal.distCal()
    except ValueError as err:
       print '> WARNING: %s is NOT calibrated: %s ' %(fitsImage,err)
if calibrate:
    #just choose a random file in the original folder in order to call the function
    dummyList = os.listdir(fdir)
    print dummyList
    firstDummy = dummyList[0]
    cal= StarCalibration(firstDummy,tfitsTable,tfitsImage,manual,paramFile=None,caldir=caldir,fdir=fdir,sedir=sedir,manCat=catOption,manCatFile=manCatFile)
    cal.distCal(addFiles=dummyList[1:])

'''
#testing scripts
#convert world coordinate(in degrees) to ARCONS coordinate
worldCoor = [98.172398,-0.0315900]
#worldCoor = [98.169492,-0.03306112]
#guide stars 20121207/112636.fits
worldCoor = [104.95365,14.241674]
worldCoor = [104.9578,14.241021]
photon = [35.9084,32.5359]
test = radec(tolError=1000)
nlist = test.centroid(worldCoor=worldCoor)
mapp = test.photonMapping('090001',15.72,14.65)
'''