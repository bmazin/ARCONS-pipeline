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
pos = 'crab nebula'
#RA,DEC = convert([['6:59:49.587','+14:14:02.96']])
#pos = RA+','+DEC

#source of data:'USNO-B1.0' or '2MASS' are usually enough. For full list: http://cdsarc.u-strasbg.fr/viz-bin/vizHelp?cats/U.htx
#source = 'USNO-B1.0'
source = '2MASS'

#name of the saved files
tfitsTable = 'test.fits'
tfitsImage = 'test_image.fits'

#it will overwrite any existing files with the same names, the file 'test_vo.xml' is not important and can be ignored
queryVizier(tfitsTable,source=source,pos=pos)
queryFitsImage(tfitsImage,'test_vo.xml',pos=pos)

#saving directory of all the calibrated files in relative path
caldir = './cal/'
#directory of fits images to be calibrated, put all the files here
fdir = './origin/'
sedir = './config/'

#if manual = False, the program will use sextractor to find source and match the correponding stars in the images
#also make sure the ./origin/ folder has appropriate sextractor parameters files and parameters
manual = False

#next, if automatic calibration is chosen, it is best to first manually correct the reference pixel coordinate on the header. This greatly increases the chances of calibrating.
refFix = False
if refFix:
    #specify the correct reference pixels. This can be found by manually open the fits file and changing the header keywords CRPIX1 AND CRPIX2 until the stars matched the one found in the catalog. This only has be done once per calibration since we assume the offset is constant.
    CRPIX1 = 560
    CRPIX2 = 360


#perform linear and polynomial calibration to each file in dir specified
for fitsImage in os.listdir(fdir):
    #fix reference value if refFix is True
    if refFix:
        #open the fits image and change CRPIX value. Save the file to *fixed.fits in the fdir(origin in defualt)
        imageList = pyfits.open(fitsImageName)
        header = imageList[0].header
        header.update('CRPIX1',CRPIX1)
        header.update('CRPIX2',CRPIX2)
        fitsImage = '%sfixed.fits' %(fitsImageName[0:6])
        try:
            os.remove(fitsImage)
        except:
            pass
        imageList.writeto(fitsImage,output_verify='ignore')        
	cal = StarCalibration(fitsImage,tfitsTable,tfitsImage,manual,paramFile=None,caldir=caldir,fdir=fdir,sedir=sedir)
	cal.linCal()
	#cal.distCal()

#convert world coordinate(in degrees) to ARCONS coordinate
worldCoor = [98.172422,-0.031578365]
test = radec(tolError=10000)
nlist = test.centroid(worldCoor=worldCoor)
