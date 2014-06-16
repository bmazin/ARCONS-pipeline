#!/bin/python
'''
Author: Yin Lin        Date: September 23, 2013

The class object StarClass loads stars from the catalog.py and from the source extractor cat file and converts them into proper format.
The lists of stars are then passed to the subclass, StarCalibration, to cross match two lists and perform offset and rotational calibration 
by changing the header keywords of the corresponding arguments. In order to use the code, you only have to call StarCalibration which inherits 
from StarClass.

'''

from functions import *
import pyfits
import matplotlib.pyplot as plt
import numpy as np
from astropy import wcs
from scipy.optimize import leastsq,fsolve
import warnings
import os
import subprocess
from glob import *
import time
from ds9 import *
import PyGuide as pg

#ignore the warning caused by astropy
warnings.filterwarnings("ignore")
    
class StarClass(object):
    """Load a list of star in the fits image and initilize stars from catalog if specified"""

    def __init__(self,fitsImageName,fitsTableName=None,manual=False,fitsCatalogName=None,caldir='./cal/',fdir='./origin/',sedir='./config/',manCat=None,manCatFile=None):
        """
        Keyword arguments:
        fitsImageName -- the input fits image file
        fitsTableName -- the input fits table file from catalog.py. Note that it has to be in the default format(asu-fits)
        manual -- if False, the program uses Sextractor to extract stars in the image, else, the program uses ds9 and PyGuide 
                  to manually calibrate the image. Also, if True, fitsCatalogName must be specified.
        fitsCatalogName -- the catalog image from catalog for manual calibration
        """

        self.fitsImageName = fitsImageName
        self.fitsTableName = fitsTableName
        self.fitsCatalogName = fitsCatalogName
        self.catalog = []
        self.starList = []
        self.sedir = sedir
        self.fdir = fdir
        self.caldir = caldir
        self.fitsdir = self.fdir + self.fitsImageName
        self.manual = manual
        self.minError = 2000
        self.calibrate = True
        self.manCat = manCat
        self.manCatFile = manCatFile
        
        #test to see if calibration is necessary
        imageList = pyfits.open(self.fitsdir)
        header = imageList[0].header
        
        try:
            #if calibrated file with suffix _offCal_rotCal already exists, skip the calibration 
            imageListCal = pyfits.open(self.caldir+self.fitsImageName[:-5]+'_offCal_rotCal.fits')
            headerCal = imageListCal[0].header
            CALERR = headerCal['CALERR']
            print 'CALERR:' + str(CALERR)  
            if CALERR > self.minError:
                raise
            self.calibrate = False
        except:
            pass
        
        if self.manual and self.calibrate :
            #defining basic parameters for PyGuide.findStars. More info can be found on help(pyGuide.findStars)
            ccd = pg.CCDInfo(0,0.00001,1,2500)

            satMask = np.zeros(image.shape)
            mask = np.zeros(image.shape)

            #this returns Centroid class instance in PyGuide
            centroidData = pg.findStars(image,mask=mask,satMask=satMask,rad=30,ccdInfo=ccd,doDS9=False)[0]

            #sort stars and discard the one that does not make sense, e.g, the stars that are out of array
            starList = []
            _count = 0
            for data in range(len(centroidData)):    
                if centroidData[_count].xyCtr[0] > 0 and centroidData[_count].xyCtr[1] > 0:
                    starList.append(centroidData[_count].xyCtr)
                _count += 1    
            
        if not self.manual and self.calibrate:
            
            #use source extractor to extractor sources, visit man page for more info
            catName = self.caldir + self.fitsImageName[:-5] + '.cat'
            paramdir = self.sedir + 'default.sex'
            checkimg = self.caldir + self.fitsImageName[:-5] + '.check'
            proc = subprocess.Popen(['sex',self.fdir+fitsImageName,'-c',paramdir,'-CATALOG_NAME',catName,'-CHECKIMAGE_NAME',checkimg])
            proc.communicate()
            #read cat file to extract stars
            starList = readCat(catName)
        
        #coordinates in pixels
        if self.calibrate:
            self.starList = np.array(starList)
    
        #if fits table is provided, loading fits table and convert the list of star into standard numpy array
        if self.manCat == 'full':
            catalog = []
            catalogAppend = readCat(self.manCatFile)
            for star in catalogeAppend:
                catalog.append([star[0],star[1]])
        
        elif self.fitsTableName != None:
            tableList = pyfits.open(self.fitsTableName)
            
            h1 = tableList[1]
            catalog = h1.data
            catalog_x = []
            catalog_y = []
            for i in range(len(catalog)):
                catalog_x.append(catalog[i][0])
                catalog_y.append(catalog[i][1])
            
            #if semi manual catalog is used, add stars in manCatFile into the existing catalog
            if self.manCat == 'semi':
                for star in readCat(self.manCatFile):
                    catalog_x.append(star[0])
                    catalog_y.append(star[1])
            catalog = zip(catalog_x,catalog_y)
        
        self.catalog = np.array(catalog)
        print len(self.catalog)
                
        #always remember to close the file to avoid memory leakage
        imageList.close()
        
    '''
    def showCatalogPlot(self):
        #show a plot of stars in catalog
        
        if self.fitsTableName == None:
        
            print 'No Fits Table Provided!'
            
        else:       
            x = []
            y = []
            for i in range(len(self.catalog)):
                #print catalog[i][0],catalog[i][1]
                x.append(self.catalog[i][0])
                y.append(self.catalog[i][1])
            plt.plot(x,y,'ro')
            plt.show()      
            
    def showImagePlot(self):
        #show a plog of stars in fits image found by the pyGuide

        starx = []
        stary = []
        for i in range(len(self.starList)):
            starx.append(self.starList[i][0])
            stary.append(self.starList[i][1])
        plt.plot(starx,stary,'ro')
        plt.gca().invert_xaxis()    
        plt.show()
    '''

class StarCalibration(StarClass):
    #Calibration utility for fits image. Catalog data needs to be provided

    def __init__(self,fitsImageName,fitsTableName,fitsCatalogName=None,manual=False,paramFile=None,caldir='./cal/',fdir='./origin/',sedir='./config/',height=3,manCat=None,manCatFile=None):
        
        #Initialize class attributes
        self.paramFile = paramFile 
        self.calibratedStar = []
        self.fitsImageName = fitsImageName
        self.fitsTableName = fitsTableName
        self.fitsCatalogName = fitsCatalogName
        self.calpix = []
        self.pix = []
        self.labelOn = True
        self.sedir = sedir
        self.caldir = caldir
        self.fdir = fdir
        self.fitsdir = self.fdir + self.fitsImageName
        self.height = height
        self.manual = manual
        
		#create folder for calibrated files if doesnt exit
        if not os.path.exists(self.caldir):
            os.makedirs(self.caldir)		
		
		#ignore the warning caused by astropy
        warnings.filterwarnings("ignore")
        
        #Initialize the super class, StarClass
        super(StarCalibration,self).__init__(self.fitsImageName,self.fitsTableName,self.manual,self.fitsCatalogName,self.caldir,self.fdir,self.sedir,manCat=manCat,manCatFile=manCatFile)
        
        #calibration needs at least 2 star to perform. Raise error if it has fewer than 2 star detected. This becomes a warning in manage.py
        if len(self.starList) < 2 and self.calibrate:        
            raise ValueError, '2 of more stars required to calculate calibration parameters. Only %s star/s detected!' %len(self.starList)
        
        #Initialize reference pixels from the header
        imageList = pyfits.open(self.fitsdir)
        header = imageList[0].header
        x1ref = header['CRPIX1']
        x2ref = header['CRPIX2']
        imageList.close()
        
        #load the star from header if it is well calibrated already
        if not self.calibrate:
            imageListCal = pyfits.open(self.caldir+self.fitsImageName[:-5]+'_offCal_rotCal.fits')
            headerCal = imageListCal[0].header
            pix = []
            calWorld = []
            calX = map(float,headerCal['CALX'].split(','))
            calY = map(float,headerCal['CALY'].split(','))
            pixX = map(float,headerCal['STARX'].split(','))
            pixY = map(float,headerCal['STARY'].split(','))         
            for count in range(len(calX)):
                pix.append([pixX[count],pixY[count]])
                calWorld.append([calX[count],calY[count]])
            self.pix = np.array(pix)
            self.calibratedStar = np.array(calWorld)
            imageListCal.close()
            
        
        #Manual cross match
        if self.calibrate and self.manual == True:
                   
            #try to see any calibrated file so no need to repeat calibration    
            try:    
                calDict = np.load(self.caldir+fitsImageName[:-5]+'calList.npz')
                self.calibratedStar = calDict['calibratedStar']
                self.pix = calDict['starListPix']         
                #Prevent python errors of too many files opened
                calDict.close()
                
            except:          
                #convert starList into world coordinate              
                self.starList = pix2world(self.fitsdir,self.starList)
                openCatalog = '-fits' + ' ' + self.fitsCatalogName

                #initialize the centroid of set of stars in catalog table and put an 'X' mark in that position
                cenX = [p[0] for p in self.catalog]
                cenY = [p[1] for p in self.catalog]
                self.centroid = [[float(sum(cenX)) / len(self.catalog), float(sum(cenY)) / len(self.catalog)]]
                
                #convert lists to standard format in order to send commands to ds9
                self.starList = convert(self.starList)
                self.catalog = convert(self.catalog)
                self.centroid = convert(self.centroid)
                
                #open ds9 instance for catalog image
                catalog = ds9(target='catalog',start=openCatalog)
                catalog.set('scale mode 99')
        
                def _turnOnLabel():
                    starCounter = 0
                    for stars in range(len(self.catalog)):
                        RA = self.catalog[stars][0]
                        DEC = self.catalog[stars][1]
                        
                        #setting the parameters for the drawings 
                        radius = 2
                        coor = 'image;'+ ' ' + 'circle' + ' ' + RA + ' ' + DEC + ' ' + str(radius)
                        catalog.set('regions',coor)
                        text = 'image; text %s %s #text="%s" font="times 15 bold"' %(RA,DEC,starCounter)
                        catalog.set('regions',text )
                        starCounter += 1
                    starCounter -= 1
                    return starCounter
            
                starCounter = _turnOnLabel()
        
                text = 'image; text %s %s #text="%s" font="times 15 bold"' %(self.centroid[0][0],self.centroid[0][1],'X')
                catalog.set('regions',text )
                
                #open ds9 instance for fits image
                openGuide = '-fits' + ' ' + self.fitsdir
                guideImage = ds9(target='guidImage',start=openGuide)
                guideImage.set('scale mode 99')
                guideImage.set('cmap value 0.9 0.7')
                guideImage.set('zoom to fit')
            
                #to account for deleting the star(key=-1), if repeat is True, the index stays the same as we move all the items 1 index to the left when we delete an item
                global repeat
                repeat = False
                #keep a list of assigned stars
                _starNumber = []
            
                count = 0
                for stars in range(len(self.starList)):
                
                    if repeat:
                    #determine whether we have to repeat the same count/index or not(due to a star being deleted by key=-1)
                        count -= 1
                        repeat = False
                
                    print count
                
                    #this try statement is to handle the special case in which we use key=-1 to delete the last star in the self.starList
                    try:
                        RA = self.starList[count][0]    
                        DEC = self.starList[count][1]
                    except:
                        break
                    
                    radius = 10
                    coor = 'image;'+ ' ' + 'circle' + ' ' + RA + ' ' + DEC + ' ' + str(radius)
                    guideImage.set('regions',coor)
    
                    while True:
    
                        def _calibration():
            
                            key = raw_input('--> ')
                
                            #the user needs to manually identify the star from the image to the one in the catalog and enter the star number assigned on the catalog. If key=-1, the user is unable to identity and star and the porgram will delete it from the list; if key=-2, the star labeling will be turning off/on, leaving only the centroid mark 'X'.     
                            try:
                                key = int(key)
                            except:
                                print 'not a number!'
                                return True
            
                            try:
                                if (key <= starCounter and key >= 0) or key == -1 or key == -2:
                                    pass    
                                else:
                                    raise ValueError
                            except:
                                print 'not within range!'
                                return True
            
                            try:
                                for no in _starNumber:
                                    if key == no:
                                        raise ValueError
                            except:
                                print 'star already assigned!'
                                return True
            
                            try:
                                if key == -1:
                                    raise ValueError
                            except:    
                                print 'unable to locate the star, skip to the next one'
                                del self.starList[count]
                                #global statement in order to make repeat variable visible outside of the function _calibration.
                                global repeat
                                repeat = True
                                return False
                        
                            try:
                                if key == -2:
                                    raise ValueError
                            except:
                                if self.labelOn:
                                    catalog.set('regions delete all')
                            
                                    #delete every label except the center 
                                    text = 'image; text %s %s #text="%s" font="times 15 bold"' %(self.centroid[0][0],self.centroid[0][1],'X')
                                    catalog.set('regions',text )
                                    self.labelOn = False
                                    return True
                                else:
                                    _turnOnLabel()
                                    self.labelOn = True
                                    return True
            
                            RA = self.catalog[key][0]
                            DEC = self.catalog[key][1]            
                            self.calibratedStar.append([RA,DEC])
                            print RA,DEC
                            _starNumber.append(key) 
            
                            #guideImage.set('regions delete all')
                            return False
        
                        if not _calibration():
                            count += 1
                            break
            
                #the loop allows for manually selecting stars on the image and specify its corresponding star number    
                while True:
                    yn = raw_input('manual selection y/n? ')
                    if yn == 'y':
                        RA = raw_input('--> RA ')
                        DEC = raw_input('--> DEC ')
                        self.starList.append([RA,DEC])
        
                        no = raw_input('--> starnumber ')
                        RA = self.catalog[int(no)][0]
                        DEC = self.catalog[int(no)][1]            
                        self.calibratedStar.append([RA,DEC])
                        print RA,DEC    
                        _starNumber.append(int(no))
                    elif yn == 'n':
                        break
                    else:
                        print 'wrong key, try again'
                
                guideImage.set('exit')
                catalog.set('exit')
        
        
                #convert from standard to degree in order to pass on to world2pix later on
                self.calibratedStar = np.array(convert(self.calibratedStar))
                self.starList = np.array(convert(self.starList))
                self.pix = world2pix(self.fitsdir,self.starList)
                
                minError = linCal()
                
                '''
                #save the calibrated list for class method reference
                saveName = self.caldir + self.fitsImageName[:-5] + 'calList'        
                np.savez(saveName,calibratedStar=self.calibratedStar,starListPix=self.pix)
                '''
                
        #automatic cross match
        if self.calibrate and self.manual == False:
            
            #convert cataglogue world coorodinates into pixel coordinates
            self.catalog = world2pix(self.fitsdir,self.catalog,self.paramFile,x1ref,x2ref)
            
            '''
            #Test to see if paramFile is provided, whether the reverse distortion transformation will give back itself as distortion is zero at reference point
            print 'test ref1:',x1ref,x2ref           
            testList = pix2world(self.fitsdir,np.array([[x1ref,x2ref]]),self.paramFile,x1ref,x2ref)                  
            print 'test ref2:', world2pix(self.fitsdir,testList,self.paramFile,x1ref,x2ref)
            raise ValueError, 'test terminated'
            '''
        
            def _patternGeneration(index,length,height): 
                """
                Generate pattern for cross matching. 
                Example 1, index=5, length=3, height=5 will return [1,1,5].
                Example 2, index=3, length=4, height=2 will return [1,1,2,1]
                Example 3, index=10, length=3, height=3 will return [3,1,1]
                
                Note that index has to start from 1.
                """
                
                height = int(height)
                index = int(index)
                length = int(length)
                
                if index > height**length:
                    raise ValueError('Index out of range!')
                else:
                    coeffList = []
                    remainder = index
                    for no in list(reversed(range(length))):
                        if no == 0:
                            coeff = remainder
                            coeffList.append(coeff)
                            break
                        coeff = np.ceil(remainder / float(height**no)) 
                        remainder = remainder - (coeff - 1)*(height**no)
                        coeffList.append(coeff)
                return coeffList
                

            #now we have both self.catalog and self.starList both in pixel coordinate, we want to cross match two lists by looking for least squares
            height = self.height        
            minList = []
            print 'total number of stars = %s' %(len(self.starList))
            
            #funciton that calculates the distance between two coordinates.
            d = lambda c1,c2:np.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2)
            
            #for each star in starList, we calculate the distance between the star and any other stars in the catalog and sort them in ascending distance in sortMinSubList
            for star in self.starList:
                minSubList = []
                index = 0
                for refStar in self.catalog:
                    minSubList.append([index,d(refStar,star)])
                    index += 1

                sortMinSubList = []
                #sort the list in ascending distance difference
                while len(sortMinSubList) < height:
                    minDist = 1000000000000
                    minIndex = 0
                    delIndex = 0
                    counter = 0
                    for param in minSubList:
                        if param[1] < minDist:
                            minIndex = param[0]
                            minDist = param[1]
                            delIndex = counter
                        counter += 1
                    sortMinSubList.append([minIndex,minDist])
                    
                    del minSubList[delIndex]
                    
                    
                    #this was meant to eliminate the identification with unreasonably far stars     
                    tolerancePixel = 150
                    count = 0
                    
                    for item in sortMinSubList:
                        #set any entries that have distance differences greater than the tolerance to None and ignore them. and not count==0 statement prevent the situation where all the entries are None.
                        if item[1] > tolerancePixel and not count == 0:
                            item[1] = None
                        count = count + 1     
                        
                minList.append(sortMinSubList)
                    
            def _matchPattern(pattern):
                sortCatalog = []         
                index = 0 
                con = False
                
                #make the starting point to be 0 compared to 1 as in _patternGeneration function
                for order in pattern:
                    order = int(order - 1)
                    #ignore None entry which has distance difference greater than the tolerance
                    if minList[index][order][1] == None:
                        sortCatalog.append(None)
                        con = True
                    else:
                        sortCatalog.append(self.catalog[minList[index][order][0]])                      
                    index += 1
                return sortCatalog,con
                    
            #just to initialize variables
            minError = 10**10
            minIndex = 1
            #loop through all possible patterns to determine which gives the least error by least squares fit        
            for index in range(height**(len(self.starList))):
                index = index + 1
                tempCatalog = self.catalog
                pattern = _patternGeneration(index,len(self.starList),height)
                sortCatalog,con = _matchPattern(pattern)
                
                #if None entry exists, skip straight to the next iteration in for loop
                if con:
                    continue
            
                self.calibratedStar = pix2world(self.fitsdir,np.array(sortCatalog),self.paramFile,x1ref,x2ref)            
                self.pix = self.starList
                
                self._offCal()
                self._rotCal()
                for i in range(1):
                    self._offCal(CD=False,openName=self.caldir + self.fitsImageName[:-5] + '_offCal_rotCal.fits')
                    error = self._rotCal(openName=self.caldir + self.fitsImageName[:-5] + '_offCal.fits')
                
                if error < minError:
                    minError = error
                    minIndex = index
     
            minPattern = _patternGeneration(minIndex,len(self.starList),height)
            print 'minimum Error is %s' %minError
            print 'with pattern %s' %minPattern
            sortCatalog = _matchPattern(minPattern)[0]
            
            #give the ordered pair of catalog stars in world coordinates(degrees) and image star in pixel coordinate which then can be passed to offCal and rotCal methods. The position of one list matches the other.
            self.calibratedStar = pix2world(self.fitsdir,np.array(sortCatalog),self.paramFile,x1ref,x2ref)            
            self.pix = np.array(self.starList)  
                        
        
        
        
        
       
       
        #create header entries to record calibrated star and error. Only perform this when the calibration is done on the first time
        if self.calibrate:
            appendCALX = []
            appendCALY = []
            appendX = []
            appendY = []
            for pixCoor in self.pix:
                appendX.append(str(pixCoor[0]))
                appendY.append(str(pixCoor[1]))
            updateHeader(self.fitsdir,'STARX',",".join(appendX))
            updateHeader(self.fitsdir,'STARY',",".join(appendY))       
            for calDeg in self.calibratedStar:
                appendCALX.append(str(calDeg[0]))
                appendCALY.append(str(calDeg[1]))
            updateHeader(self.fitsdir,'CALX',",".join(appendCALX))
            updateHeader(self.fitsdir,'CALY',",".join(appendCALY))
            updateHeader(self.fitsdir,'CALERR',minError)
            
             
        
    def linCal(self,iteration=15):
        """
        A wrapper around offCal and rotCal methods. It performs translational and rotational calibration in appropriate order. The calibration will be performed until either tolerance error or iteration upper bound is met. Iteration argument is obsolete. 
        
        Keyword arguments:
        NONE, iteration argument is obsolete.
        """
        
        if not self.calibrate:
            imageList = pyfits.open(self.caldir+self.fitsImageName[:-5]+'_offCal_rotCal.fits')
            header = imageList[0].header
            CALERR = header['CALERR']
            return CALERR
        
        self._offCal()
        error = self._rotCal()
        upperBound = 30
        tolerance = 1
        errorTemp = 1000000000
        i = 0
        
        while i < 30:
            self._offCal(CD=False,openName=self.caldir + self.fitsImageName[:-5] + '_offCal_rotCal.fits')
            error = self._rotCal(openName=self.caldir + self.fitsImageName[:-5] + '_offCal.fits')
            if abs(errorTemp - error) < tolerance:
                break
            errorTemp = error
            i += 1
        print 'calibration finished in %s iterations!' %i
            
        print 'error = %s' %error
        return error
        '''
        #record the name and error in runRecord.txt
        recordFile = 'runRecord.txt'      

        try:
            os.remove(recordFile)
        except: 
            pass

        openFile = open(recordFile,'a')
        outputList = ''.join([self.fitsImageName.ljust(30),str(error/len(self.starList)).ljust(30),'\n'])
        openFile.write(outputList)
        openFile.close()  
        '''  
                
    def _offCal(self,CD=True,openName=None):
        """
        Perform translational offset calibration by calibrating the reference pixel WCS coordinate
        
        Keyword arugments:
        CD --- if True, the program will recompute the CD matrix based on the catalog and vice versa (default:True)
        openName --- the name of the fits file to be calibrated (default:self.fitsImageName[:-5]+'offCal.fits'). This must match the file of self.fitsImageName.  
        """
       
        #default arguement
        if openName == None:
            openName = self.fitsdir

        imageList = pyfits.open(openName)
        header = imageList[0].header
        if header['CTYPE1'] == 'RA--TAN':
            header['CTYPE1'] = 'RA---TAN'
        
        x1ref = header['CRPIX1']
        x2ref = header['CRPIX2']
        polydeg = 'dist'
        
        #calculate new CD matrix from calibrated star,this should be initilize only for the first offCal,the negative of the first componenet comes from the fact the x-axis(East direction) is inverted in our images
        if CD:       
            header['CD1_1'] = -header['PLTSCALE']
            header['CD2_2'] = header['PLTSCALE'] 
            header['CD1_2'] = 0
            header['CD2_1'] = 0
            
        #each time when we call different calibration class methods, we have to convert the wcs of catalog star into pix using DIFFERENT image headers which causes shift in pixel coordinates!!
        self.calpix = world2pix(openName,self.calibratedStar,self.paramFile,x1ref,x2ref)
        
        #just to make sure they are in numpy array so we can pass it to functions below
        pix = np.array(self.pix)
        calpix = np.array(self.calpix)
           
        '''
        def poly2(p,x1,x2):
            a00,a01,a02,a10,a11,a12,a20,a21,a22 = p
            y = a00+a01*x2+a02*x2**2+a10*x1+a11*x1*x2+a12*x1*x2**2+a20*x1**2+a21*x1**2*x2+a22*x1**2*x2**2
            return y
            
        def poly1(p,x1,x2):
            a00,a01,a10,a11,a20,a02 = p
            y = a00 + a01*x2 + a10*x1 + a11*x1*x2 + a20*x1**2 + a02*x2**2
            return y
   
        def leastDist(p,x1,x2):
            x0,y0 = p
            y = np.sqrt((x1+x0)**2+(x2+y0)**2)

     
        def residuals(p,y,x1,x2,polydeg=2):
            if polydeg == 2:
                err = y - poly2(p,x1,x2)
            elif polydeg == 1:
                err = y - poly1(p,x1,x2)
            return err
            
        '''
        #residual of least squares    
        def residuals(p,y1,y2,x1,x2):
            x0,y0 = p
            err = np.sqrt((y1-(x1+x0))**2 + (y2-(x2+y0))**2)
            return err
                
        #initialize starting parameters, doesn't matter where it starts 
        if polydeg == 2:
            p0 = [0.001]*9
        elif polydeg == 1:
            p0 = [0.001]*6
        elif polydeg == 'dist':
            p0 = [100]*2
        
        '''
        starDist = []
        for star in self.pix:
            dist = (-header['CRPIX1'] + star[0])**2 + (star[1] - 844 + header['CRPIX2'])**2     
            starDist.append(dist)
        #[minimum position,element]
        mini = [0,starDist[0]]
        for count in range(1,len(starDist)):
            if starDist[count] < mini[1]:
                mini = [count,starDist[count]]      
        header.update('CRPIX1',self.pix[mini[0]][0])
        header.update('CRPIX2',self.pix[mini[0]][1])
    
        x1ref = header['CRPIX1']
        x2ref = header['CRPIX2']
        '''
        
        x1 = np.array([(no) for no in (pix[cr][0] for cr in range(len(pix)))])
        x2 = np.array([(no) for no in (pix[cr][1] for cr in range(len(pix)))])
        '''
        if not polydeg == 'dist':
            dtor = calpix - pix
            y1meas = np.array([no for no in (dtor[cr][0] for cr in range(len(dtor)))])
            y2meas = np.array([no for no in (dtor[cr][1] for cr in range(len(dtor)))])  
            dt1 = leastsq(residuals,p0,args=(y1meas,x1,x2,polydeg))
            dt2 = leastsq(residuals,p0,args=(y2meas,x1,x2,polydeg))
            np.savez('polyparams',dt1=dt1,dt2=dt2)
        
            #calculate shift in wcs of reference cooridnate and apply to header
            x1ref = x1ref + dt1[0][0]
            x2ref = x2ref + dt2[0][0]
        '''
        
        #least squares fits to minimize distance difference between catalog and our image
        y1 = np.array([(no) for no in (calpix[cr][0] for cr in range(len(calpix)))])
        y2 = np.array([(no) for no in (calpix[cr][1] for cr in range(len(calpix)))])
        
        #With fulloutput=True, we can obtain the residual/error of the least squares fit
        dt = leastsq(residuals,p0,args=(y1,y2,x1,x2))
        
        #update the header
        x1refNew = x1ref + dt[0][0]
        x2refNew = x2ref + dt[0][1] 
        refPix = [[x1refNew,x2refNew]]
        
        #calculate appropriate CRVALS and update header
        world = pix2world(openName,refPix,paramFile=self.paramFile,crval1=x1ref,crval2=x2ref)                  
        header['CRVAL1'] = world[0][0]
        header['CRVAL2'] = world[0][1] 
               
        #remove existing file and save the calibrated file 
        saveName = self.caldir + self.fitsImageName[:-5] + '_offCal.fits'
        try:
            os.remove(saveName)
        except:
            pass
        imageList.writeto(saveName,output_verify='ignore')
                
    
    def _rotCal(self,openName=None):
        """
        perform rotational calibration of fits image by multiplying a rotational matrix onto the CD matrix which then creates a new CD matrix.
        
        Keyword arguments:
        openName --- the name of the fits file to be calibrated (default:self.fitsImageName[:-5]+'offCal.fits'). This must match the file of self.fitsImageName.      
        """
        
        if openName == None:    
            openName = self.caldir + self.fitsImageName[:-5] + '_offCal.fits'
            
        #apply rotCal after offCal!
        imageList = pyfits.open(openName)
        header = imageList[0].header
            
        x1ref = header['CRPIX1']
        x2ref = header['CRPIX2']
        
        self.calpix = world2pix(openName,self.calibratedStar,self.paramFile,x1ref,x2ref)
        
        pix = np.array(self.pix)
        calpix = np.array(self.calpix)
        
        #initialize CD matrix
        CD11 = header['CD1_1']
        CD12 = header['CD1_2']
        CD21 = header['CD2_1']
        CD22 = header['CD2_2']
        
        CD = np.matrix([[CD11,CD12],[CD21,CD22]])

        x1ref = header['CRPIX1']
        x2ref = header['CRPIX2'] 
        xref = np.array([x1ref,x2ref])
        
        #first component of vector after being rotated
        def y1m(p,x1,x2):           
            theta = p           
            y = np.cos(theta)*x1 - np.sin(theta)*x2
            return y
            
        #second component of vector after being rotated
        def y2m(p,x1,x2):
            theta = p
            y = np.sin(theta)*x1 + np.cos(theta)*x2
            return y

        #residuals of least squares fit
        def residuals(p,y1,y2,x1,x2):
            err = np.sqrt((y1 - y1m(p,x1,x2))**2 + (y2-y2m(p,x1,x2))**2)
            return err

        #initialize guessing parameters, doesn't matter
        p0 = [0.1]
        
        '''
        #I think this line is acutally false, let me try
        #if paramFile is provided, we need to first apply distortion prior to rotCal()
        if self.paramFile != None:
            pix = distApp(pix,self.paramFile,x1ref,x2ref)
        '''
        
        x1 = np.array([(no-x1ref) for no in (pix[cr][0] for cr in range(len(pix)))])
        x2 = np.array([(no-x2ref) for no in (pix[cr][1] for cr in range(len(pix)))])
        y1meas = np.array([(cal-x1ref) for cal in (calpix[cr][0] for cr in range(len(calpix)))])
        y2meas = np.array([(cal-x2ref) for cal in (calpix[cr][1] for cr in range(len(calpix)))])
        
        #With fulloutput=True, we can obtain the residual/error of the least squares fit
        dt,junk1,infoDict,junk2,junk3 = leastsq(residuals,p0,args=(y1meas,y2meas,x1,x2),full_output=True)
        error = (infoDict['fvec']**2).sum()
                
        #rotational matrix components (counter-clockwise)
        theta = dt[0]
        R11 = np.cos(theta)
        R12 = -np.sin(theta)
        R21 = np.sin(theta)
        R22 = np.cos(theta)
        R = np.matrix([[R11,R12],[R21,R22]])
        
        #print 'theta = %s' %(theta)
        
        #multiply CD matrix by the rotation matrix, R. Remember we apply CD first then we apply R(or the other way depend on how you compute R)!
        CD = np.dot(CD,R)
        
        #turn matrix into a list to extract components
        CD = CD.tolist()     
        
        #updating header keywords  
        for p in range(2):
            for q in range(2):
                keyword = 'CD%s_%s' %(p+1,q+1)
                header[keyword] = CD[p][q]
            
        #remove existing file and save the calibrated image
        saveName = self.caldir + self.fitsImageName[:-5] + '_offCal_rotCal.fits'       
        try:
            os.remove(saveName)
        except:
            pass
        
        imageList.writeto(saveName,output_verify='ignore') 
        
        #return the total residual of least squares fit to know the quality of fitting  
        return error
        
    def distCal(self,openName=None,addFiles=[]):
        """
        Perform distortion calibration using SIP convention (http://fits.gsfc.nasa.gov/registry/sip.html) on the fits image. If paramFile argument is provided, it simply applies the parameters
        from the param file without doing actual calibration.
        
        Keyword arugments:
        openName --- the name of the fits file to be calibrated (default:self.fitsImageName[:-5]+'offCal.fits'). This must match the file of self.fitsImageName
        addFiles --- any additional fits image that will be calibrated with openName to provide more data point. Note that the images provided here have to be calibrated priorly (default:None)
        """        
        #Set default arguments
        if openName == None:    
            openName = self.caldir + self.fitsImageName[:-5] + '_offCal_rotCal.fits'
        
        #Initialize reference pixels from the header
        imageList = pyfits.open(openName)
        header = imageList[0].header
        x1ref = header['CRPIX1']
        x2ref = header['CRPIX2']
            
        #Calculate distortion parametes if paramFile is missing
        if self.paramFile == None:       
            self.calpix = world2pix(openName,self.calibratedStar)
            if addFiles != None:
                #add the directory of add files
                for tfile in addFiles:
                    #loading calibrated star from header
                    addImageList = pyfits.open(self.fdir+tfile)
                    addheader = addImageList[0].header
                    pix = []
                    calWorld = []
                    
                    #append the pix and corresponding stars into a single list. Then, the calibration parameters will base on these stars.
                    calX = map(float,addheader['CALX'].split(','))
                    calY = map(float,addheader['CALY'].split(','))
                    pixX = map(float,addheader['STARX'].split(','))
                    pixY = map(float,addheader['STARY'].split(','))         
                    for count in range(len(calX)):
                        pix.append([pixX[count],pixY[count]])
                        calWorld.append([calX[count],calY[count]])
                    appendPix = np.array(pix)
                    appendCal = world2pix(self.caldir+tfile[:-5]+'_offCal_rotCal.fits',np.array(calWorld))
                    '''
                    if self.manual:
                        nlist = np.load(tfile[:-5]+'calList.npz')
                        appendPix = nlist['starListPix']
                    elif not self.manual:
                        catName = tfile[:-5] + '.cat'
                        appendPix = readCat(catName)
                        #problems
                    
                    appendCal = world2pix(tfile[:-5]+'_offCal_rotCal.fits',nlist['calibratedStar']) 
                    '''  
                    try:
                        self.pix = self.pix.tolist()
                        self.calpix = self.calpix.tolist()
                    except: 
                        pass    
                    
                    for coor in appendPix.tolist():
                        self.pix.append(coor)
                    self.pix = np.array(self.pix)
                    
                    for coor in appendCal:
                        self.calpix.append(coor)
                    self.calpix = np.array(self.calpix) 
            
                    
            imageList = pyfits.open(openName)
            header = imageList[0].header
            
            #residual of least squares for using leastsq function
            def residuals(p,y,x1,x2):
                err = y - poly2(p,x1,x2)
                return err
            '''
            def polyTPV(p,x1,x2):
                PV1_0,PV1_1,PV1_2,PV1_3,PV1_4,PV1_5,PV1_6,PV1_7,PV1_8,PV1_9,PV1_10,PV1_11 = p
                r = np.sqrt(x1**2+x2**2)            
                y = PV1_0 + PV1_1*x1 + PV1_2*x2 + PV1_3 * r + PV1_4*x1**2 + PV1_5*x1*x2 + PV1_6*x2**2 + PV1_7*x1**3 + PV1_8*x1**2*x2 + PV1_9*x1*x2**2 + PV1_10*x2**3 + PV1_11*r**3
                return y
                
            def residualsTPV(p,y,x1,x2):
                err = y - polyTPV(p,x1,x2)
                return err
            '''   
            
            #initialize starting parameters, IT MATTERS WHERE IT STARTS!! CHOOSE A SMALL NUMBER FOR TPV (0.000001)  
            p0 = [0.00000001]*9
            
            x1ref = header['CRPIX1']
            x2ref = header['CRPIX2']
            xref = np.array([x1ref,x2ref])
            
            self.calpix = np.array(self.calpix)
            self.pix = np.array(self.pix)    
            
            dtor = self.calpix - self.pix
            y1meas = np.array([no for no in (dtor[cr][0] for cr in range(len(dtor)))])
            x2 = np.array([(no-x2ref) for no in (self.pix[cr][1] for cr in range(len(self.pix)))])
            y2meas = np.array([no for no in (dtor[cr][1] for cr in range(len(dtor)))])
            x1 = np.array([(no-x1ref) for no in (self.pix[cr][0] for cr in range(len(self.pix)))])

            '''         
            #initialize CD matrix
            CD11 = header['CD1_1']
            CD12 = header['CD1_2']
            CD21 = header['CD2_1']
            CD22 = header['CD2_2']
            
            CD = np.matrix([[CD11,CD12],[CD21,CD22]])
            
            interCal = []
            interPix = []
            for cal in self.calpix:
                #convert into coordinate relative to the reference coorindate since CD matrix acts on this space
                cal = np.array(cal) - xref          
                product = np.dot(CD,cal)        
                #convert a matrix class(in this case, a 2x1 vector) into a list
                product = product.tolist()[0]
                interCal.append(product)                
            for pix in self.pix:            
                pix = np.array(pix) - xref
                product = np.dot(CD,pix)        
                product = product.tolist()[0]
                interPix.append(product)
            
            x1 = np.array([no for no in (interPix[cr][0] for cr in range(len(interPix)))])
            x2 = np.array([no for no in (interPix[cr][1] for cr in range(len(interPix)))])
            y1meas = [cal for cal in (interCal[cr][0] for cr in range(len(interCal)))]
            y2meas = [cal for cal in (interCal[cr][1] for cr in range(len(interCal)))]
            '''
            print 'calculating distortion parameters...'
            dt1,junk1,infoDict,junk2,junk3 = leastsq(residuals,p0,args=(y1meas,x1,x2),full_output=True)
     
            #REMEMBER TO CHANGE THE ORDER OF X1 AND X2 WHEN SWITCH FROM SIP TO TPV!!!!!!IMPORTATNT!!!!!!
            dt2,junk1,infoDict,junk2,junk3 = leastsq(residuals,p0,args=(y2meas,x1,x2),full_output=True)
            
            #conventions used in SIP, do not change this
            coeffList = ['10','01','20','02','11','21','12','30','03']
            
            #save the parameter file (everything saved will have suffix .npz added automatically)
            psaveName = self.caldir + 'params'
            try:
                os.remove(psaveName)
            except:
                pass
            
            #This saves the distortion parameters in 'params.npz' which can be reused later
            np.savez(psaveName,dt1=dt1,dt2=dt2,xref=xref,coeffList=coeffList)
                       
            #update the header of every fits images provided
            saveName = self.caldir + self.fitsImageName[:-5]+'_allCal.fits'
            coeffUpdate(saveName,openName,dt1,dt2,coeffList)
            for tfile in addFiles:
                coeffUpdate(self.caldir+tfile[:-5]+'_allCal.fits',self.caldir+tfile[:-5]+'_offCal_rotCal.fits',dt1,dt2,coeffList)
                   
                
            '''
            header.update('CTYPE1','RA---TPV')
            header.update('CTYPE2','DEC--TPV')
            
            for no in range(len(dt1[0])):
                keyword1 = 'PV1_%s' %(no)
                keyword2 = 'PV2_%s' %(no)
                key1 = dt1[0][no]
                key2 = dt2[0][no]
                header.update(keyword1,key1)
                header.update(keyword2,key2)
            ''' 
            
        #apply coefficients to header if paramFile is provided
        if self.paramFile != None:          
            distHeaderUpdate(openName,self.fitsImageName[:-5]+'_allCal.fits',self.paramFile)
         
if __name__ == '__main__':
     
    #PSR06+14 coordinate
    #print convert([['6:59:49.587','+14:14:02.96']])

    #SDSS J0651
    #print convert([['6:51:33.338','+28:44:23.37']]) 
    
    fitsTableName = 'test.fits'
    fitsImageName = '071023fix.fits'
    fitsCatalogName = 'test_image.fits'
    paramFile = 'params.npz'
    
    cal = StarCalibration(fitsImageName,fitsTableName,fitsCatalogName,manual=True,paramFile=None)
    cal.linCal()
    cal.distCal(addFiles=['100011fix.fits'])
 
    '''   
    fileList = glob('*fix.fits')
   
    paramFile = 'runRecord.txt'
    
    open(paramFile,'w').close()

    count = 0
    startTime = time.time()
    for nfile in fileList:
     
        fitsTableName = 'test.fits'
        fitsImageName = nfile
        fitsCatalogName = 'test_image.fits'

        addFiles = []
        
    print 'solving %s...' %(nfile)
    subStartTime = time.time()
        cal = starCalibration(fitsImageName,fitsTableName,fitsCatalogName)
    subEndTime = time.time() - subStartTime
    print 'done in %s seconds' %(subEndTime)
    
    elapsedTime = time.time()-startTime
    print 'Everything finised in %s seconds' %elapsedTime
    '''
    '''
    cal.offCal(polydeg='dist',addFiles=addFiles)
    #cal.rotCal(openName=fitsImageName[:-5] + '_offCal.fits')   
    '''
    '''
    #iterative calibration of offCal and rotCal. For 071023.fits, it takes only 2 loops to converge.    
    for a in range(20):
        addFilesIt = ['crab2_offCal_rotCal.fits','crab3_offCal_rotCal.fits']
        addFilesIt = []
        cal.offCal(CD=False,openName=fitsImageName[:-5] + '_offCal_rotCal.fits',polydeg='dist',addFiles = addFilesIt)
        

        for nfile in addFiles:
            ncal = starCalibration(nfile,fitsTableName,fitsCatalogName)
            ncal.rotCal(openName=nfile[0:5]+'_offCal.fits')
        cal.rotCal(openName=fitsImageName[:-5] + '_offCal.fits')
    
    addFiles = ['crab1.fits']
    cal.distCal(addFiles=addFiles)

    
    

    '''
