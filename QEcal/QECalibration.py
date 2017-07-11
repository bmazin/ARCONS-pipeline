'''
Author: Alex Walter
Date: 5/13/2013

For computing and comparing QE measurements
(Originally for comparing QE measurements with different polariztions of light.
 If you just want to look at a single QE measurement choose a dummy polariztion angle = 0)
 
Usage:
$ python QECalibration.py
Then click the buttons

Advanced Usage:
Save the QEData.npy as something specific. ie QEData_*date*.npy
Copy the initializeWithAug22Data() function but fill in the relevant info
Call your new initialze function in the GUI __init__(). 

The way it works:
1. All the photon data is taken in a single long exposure. 
   We must first determine which section of the exposure corresponds to which wavelength of light.
   This is saved as *QEfile*_param.txt
2. Calculate the QE for each pixel. This is saved as "QEData.npy" (which can be reloaded if you restart the GUI)
3. Some pixels may be bad and skew the average, we'll flag them and ignore them. 
   Then we plot the average, median, mode QE across the array 
4. We can compare one QE measurement to another.
   For example, with different optimal filters loaded, different polarizations of light, different cool downs, etc..

Important Files:
1. obsFile.h5           - This contains the raw photon data from the QE measurement
2. QEfile.txt           - This contains the output from the QE labview code with the calibrated photodiode intensity data
3. QEfile_param.txt     - This file is created with this GUI. It contains the timestamps during the obsFile.h5 during which each wavelength was shown on the MKID array
4. QEData.npy           - This is a temporary file created by this GUI containing the calculated pixel QE data. I save this because it can take a while to do the calculation and I don't want to redo it evertime I open reopen the GUI

'''

import sys, os, time, struct, math
from os.path import isfile
from numpy import *
from tables import *
from PyQt4 import QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from mpl_pyqt4_widget import MPL_Widget
from matplotlib import pyplot as plt
import QEfile
from QEfile import *


class QECalibration(QMainWindow):
    def __init__(self, parent=None):
        QMainWindow.__init__(self, parent)
        #print '1'
        self.good_pixel = []
        self.goodPixelQELimit = 0.1
        self.selectedPixel = None
        self.obsFileList = []
        self.numObservations = 0
        self.bmapList = []              #different observations have different unix times
        self.QEFileList = []
        self.currentObs=-1
        self.setWindowTitle('QE Calibration')
        self.imagex = 308
        self.imagey = 322
        self.obsTime = 2700
        self.create_main_frame()
        self.position = self.tv_image.pos()
        #self.addtomainframe()
        self.create_status_bar()

        #use mouse to select pixel from tv_image, also triggers a new spectrum to be displayed
        self.tv_image.mousePressEvent = self.start_pixel_select

        #self.beammapfile = os.environ['BEAMMAP_PATH']#"beamimage.h5"   #Palomar 2012 device
        #self.beammapfile='sorted_beamimage46x44.h5'

        #load beam map from default beammap directory
        #self.loadbeammap()

        #self.initializeWithNov3Data()
        #self.initializeWithMay12Data()
        #self.initializeWithNov14Data()
        #self.initializeWithNov17Data()
        #self.initializeWithNov23Data()
        #self.initializeWithNov25Data()
        #self.initializeWithAug22Data()
        
    def initializeWithAug22Data(self):
        #Sci 6
        obsfiles = ['obs_20140822-184134.h5']
        QEfiles = ['QE_20140822-184134.txt']
        path = '/Scratch/QEData/20140822/'
        os.system("cp QEData_obs_20140822-184134.npy QEData.npy")

        angleList = [0]
        for i in range(len(obsfiles)):
            self.addMeasurementFile()
            #print str(path + obsfiles[i])
            self.initializeObsFile(str(path + obsfiles[i]))
            self.initializeQEFile(str(path + QEfiles[i]),int(angleList[i]))
            self.doneCal()

    def initializeWithNov25Data(self):
        #Sci 5 - optimal filters
        obsfiles = ['obs_20131125-202654.h5']
        QEfiles = ['QE_20131125-202654.txt']
#        path = '/home/abwalter/ARCONS-pipeline/QECalibration/'
        path = '/Scratch/QEData/20131125/'
        os.system("cp QEData_obs_20131125-202654.npy QEData.npy")

        angleList = [0]
        for i in range(len(obsfiles)):
            self.addMeasurementFile()
            #print str(path + obsfiles[i])
            self.initializeObsFile(str(path + obsfiles[i]))
            self.initializeQEFile(str(path + QEfiles[i]),int(angleList[i]))
            self.doneCal()

    def initializeWithNov23Data(self):
        #Sci 5
        obsfiles = ['obs_20131124-022103.h5']
        QEfiles = ['QE_20131124-022103.txt']
#        path = '/home/abwalter/ARCONS-pipeline/QECalibration/'
        path = '/Scratch/QEData/20131123/'
        os.system("cp QEData_obs_20131124-022103.npy QEData.npy")

        angleList = [0]
        for i in range(len(obsfiles)):
            self.addMeasurementFile()
            #print str(path + obsfiles[i])
            self.initializeObsFile(str(path + obsfiles[i]))
            self.initializeQEFile(str(path + QEfiles[i]),int(angleList[i]))
            self.doneCal()

    def initializeWithNov17Data(self):
        #Sci 4
        obsfiles = ['obs_20131118-025832.h5']
        QEfiles = ['obs_20131118-025832.txt']
#        path = '/home/abwalter/ARCONS-pipeline/QECalibration/'
        path = '/Scratch/QEData/20131117/'
        os.system("cp QEData_obs_20131118-025832.npy QEData.npy")        

        angleList = [0]
        for i in range(len(obsfiles)):
            self.addMeasurementFile()
            #print str(path + obsfiles[i])
            self.initializeObsFile(str(path + obsfiles[i]))
            self.initializeQEFile(str(path + QEfiles[i]),int(angleList[i]))
            self.doneCal()

    def initializeWithNov14Data(self):
        #Sci 5  (bad coax on feedline 2)
        obsfiles = ['obs_20131115-023325.h5']
        QEfiles = ['obs_20131115-023325.txt']
#        path = '/home/abwalter/ARCONS-pipeline/QECalibration/'
        path = '/Scratch/QEData/20131114/'
        os.system("cp QEData_obs_20131115-023325.npy QEData.npy")  

        angleList = [0]
        for i in range(len(obsfiles)):
            self.addMeasurementFile()
            #print str(path + obsfiles[i])
            self.initializeObsFile(str(path + obsfiles[i]))
            self.initializeQEFile(str(path + QEfiles[i]),int(angleList[i]))
            self.doneCal()

    def initializeWithNov3Data(self):
        obsfiles = ['obs_20131104-030007.h5']
        QEfiles = ['QE_20131103.txt']
#        path = '/home/abwalter/ARCONS-pipeline/QECalibration/'
        path = '/Scratch/QEData/20131103/'
        
        angleList = [0]
        for i in range(len(obsfiles)):
            self.addMeasurementFile()
            #print str(path + obsfiles[i])
            self.initializeObsFile(str(path + obsfiles[i]))
            self.initializeQEFile(str(path + QEfiles[i]),int(angleList[i]))
            self.doneCal()

    def initializeWithMay12Data(self):
        obsfiles = ['obs_20130512-215105.h5', 'obs_20130512-224055.h5','obs_20130512-232518.h5','obs_20130513-000311.h5','obs_20130513-013920.h5', 'obs_20130513-021837.h5','obs_20130513-025834.h5','obs_20130513-033611.h5','obs_20130513-041526.h5']
        QEfiles = ['20130512POL1.txt', '20130512POL2.txt', '20130512POL3.txt', '20130512POL4.txt', '20130512POL5.txt', '20130512POL6.txt', '20130512POL7.txt', '20130512POL8.txt', '20130512POL9.txt']
#        path = '/home/abwalter/ARCONS-pipeline/QECalibration/'
        path = '/Scratch/PolarizationData/'
        
        angleList = [90, 90, 45, 45, 0, 0, -45, -45, -90]
        for i in range(len(obsfiles)):
            self.addMeasurementFile()
            #print str(path + obsfiles[i])
            self.initializeObsFile(str(path + obsfiles[i]))
            self.initializeQEFile(str(path + QEfiles[i]),int(angleList[i]))
            self.doneCal()
            

    def toggleGoodPixel(self):
        if (self.selectedPixel != None) & (self.currentObs != -1):
            x = self.selectedPixel%(self.nxpix)
            y = self.selectedPixel/self.nxpix


            try:
                self.good_pixel.remove(self.selectedPixel)
                print 'Removed '+str(self.bmapList[self.currentObs][y][x])+' from list'
                self.scene.addRect(self.scalex*(x),self.scaley*((self.nypix -1)-y),(self.scalex),(self.scaley), Qt.blue)
            except ValueError:
                self.good_pixel.append(self.selectedPixel)
                print 'Added '+str(self.bmapList[self.currentObs][y][x])+' to list!'
                self.scene.addRect(self.scalex*(x),self.scaley*((self.nypix -1)-y),(self.scalex),(self.scaley),Qt.darkCyan)

    def getAvgQE(self):
        # get QE of all pixels at all wavelength for current obs
        # set bad wavelengths to -1
        # plot median
        self.getQEData()
        medianArr=[]
        avgArr=[]
        modeArr=[]
        wavelengthArr = (self.QEFileList[0]).data[:,0]
        r=0.2
        numBins=50
        
        for i in range(len(wavelengthArr)): #number of wavelengths
            data=np.copy(self.QEData[i,self.currentObs])        #44x46 array of QE's for pixels at wavelength i for current obs
            #remove roach 4,5,6,7
            #for k in range(self.nxpix):
            #    for j in range(self.nypix):
            #        if self.bmap[j,k][2:3]=='4' or self.bmap[j,k][2:3]=='5' or self.bmap[j,k][2:3]=='6' or self.bmap[j,k][2:3]=='7':
            #            data[k,j]=-1
            #ignore pixels with <0.1% QE
            #        if data[k,j]<0.001:
            #            data[k,j]=-1
            hist=np.histogram(data,bins = numBins,range = (0.0, r),density=False)
            hist[0][0]=0
            modeArr.append((np.argmax(hist[0])+0.5)*r/numBins)
            
            data=np.ma.array(data,mask=data<0.001)
            medianQE=np.ma.median(data)
            medianArr.append(medianQE)
            averageQE=np.ma.average(data)
            avgArr.append(averageQE)

        plt.plot(wavelengthArr,modeArr,label='mode')
        plt.plot(wavelengthArr,medianArr,label='median')
        plt.plot(wavelengthArr,avgArr,label='avg')
        plt.legend()
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('QE')
        plt.show()

        #np.savetxt('avgQE_20131125-202654.txt',np.transpose([wavelengthArr,medianArr]),fmt='%i\t%.4f')
            


    def compareQE(self):
        # get QE of all pixels at all wavelengths with current obs
        # set bad wavelengths to -1
        # compare other observations to this one
        # make average of each wavelength over different pixels in different obs files compared to their pixel in current obs
        # plot for different wavelengths
        # average wavelengths
        self.getQEData()

        numRows = len(self.QEData)      #number of wavelengths
        numCol = len(self.QEData[0])    #number of observation files
        plt.figure()
#        QEFileObject = self.QEFileList[self.currentObs]
        medianQEArr = np.zeros((numRows,numCol))
        angleArr = []
        wavelengthArr = (self.QEFileList[0]).data[:,0]
        for QEFileObject in self.QEFileList:
            angleArr.append(QEFileObject.angle)

        print "Making histogram plots..."
        for i in range(numRows):
            for j in range(numCol):
#                plt.subplot(numRows,numCol,i*numCol+j+1)
                plt.subplot(1,numCol,0*numCol+j+1)
#                x=self.QEData[i,self.currentObs]
#                y=self.QEData[i,j]
                divided = np.select([((self.QEData[i,self.currentObs]>0) & (self.QEData[i,j]>0))],[1.0*self.QEData[i,j]/self.QEData[i,self.currentObs]],-1)
                weights = np.select([(divided>0)],[1],0)
                if np.sum(weights) != 0:                
                    averageQE = np.average(divided,None,weights)
                    medianQE = np.median(np.sort(divided.flatten())[np.where(np.sort(divided.flatten())>0)[0][0]:])
                else:
                    averageQE = -1
                    medianQE = -1
                #print '(obs/currentObs,wavelength): (' +str(j)+','+str(i)+') average relative QE: '+str(averageQE)

                #print '(obs/currentObs,wavelength): (' +str(j)+','+str(i)+') median relative QE: '+str(medianQE)
                medianQEArr[i,j] = medianQE
                plt.hist(divided.flatten(), bins = 100,range = (0.0, 2.0))
            print 'wavelength: ' + str(wavelengthArr[i])
            plt.show()

#        plt.show()

        print "Making median plots..."
        for i in range(len(wavelengthArr)):
#        for i in range(len(wavelengthArr[0:13])):
#        for i in [6]:
            plt.plot(angleArr, medianQEArr[i,:],marker='o',label=str(wavelengthArr[i]))
            plt.legend()
            plt.show()

#        plt.legend()
#        plt.show()


#        wavelength = 700
#        print "wv: " +str(wavelength)
#        x = self.QEData[0,0]
#        print x
#        print x.shape
#        y = self.QEData[0,1]
#        print y
#        print y.shape
        
#        divided = np.select([((x>0) & (y>0))],[1.0*y/x],-1)

#        weights = np.select([(divided>0)],[1],0)
#        averageQE = np.average(divided,None,weights)
#        print averageQE
        #average(self.QEData[wavelength,obs 1]/self.QEData[wavelength, obs 2])
#        plt.hist(divided.flatten(), range = (0.0, 2.0), bins = 100)
#        plt.show()

    def histogramQE(self):
        # get QE of all pixels at all wavelengths with current obs
        # set bad wavelengths to -1
        # histogram QE's for specified obs and wavelength
        self.getQEData()

        #numRows = len(self.QEData)      #number of wavelengths
        #numCol = len(self.QEData[0])    #number of observation files
        plt.figure()

        wavelengthArr = (self.QEFileList[0]).data[:,0]
        print wavelengthArr

        print 'Making histogram plots...'
        #plt.subplot(numRows,1,1)
        for i in range(len(wavelengthArr)):
            plt.subplot(np.ceil(np.sqrt(len(wavelengthArr))),np.ceil(np.sqrt(len(wavelengthArr))),i+1)
            #plt.hist(self.QEData[i,self.currentObs].flatten(),range=(0.0,0.2),bins=50)
            data=np.copy(self.QEData[i,self.currentObs])
            #remove roach 4,5,6,7
            #for k in range(self.nxpix):
            #    for j in range(self.nypix):
            #        if self.bmap[j,k][2:3]=='4' or self.bmap[j,k][2:3]=='5' or self.bmap[j,k][2:3]=='6' or self.bmap[j,k][2:3]=='7':
            #            data[k,j]=-1
            #ignore pixels with 0 QE
            #        if data[k,j]<0.0001:
            #            data[k,j]=-1
            numBad=np.sum(data<=0)
            print 'Wavelength: '+str(wavelengthArr[i])+'    numBad: '+str(numBad)
            plt.hist(data.flatten(),range=(0.0,0.2),bins=50)
            ax=plt.gca()
            #ax.legend()
            ax.set_title(str(wavelengthArr[i])+' nm')
            ax.set_ylim((0.,(self.nxpix*self.nypix - numBad)/5.))
        plt.show()

        

    def getQEData(self):
        


        if isfile('QEData.npy'):
            self.QEData = np.load('QEData.npy')
            print 'Loaded QEData from: QEData.npy'
            print "`-->Misses any recently added data!"
        else:
            if (self.currentObs>=0):           
                try:
                    self.obsFileList[self.currentObs]
                    self.QEFileList[self.currentObs]
                except IndexError:
                    print 'Missing observation or QE file'
                
                dataQE = []
                QEFileObject = self.QEFileList[self.currentObs]
                print "Populating QE Data"
                for wavelength in QEFileObject.data[:,0]:
#                for wavelength in [700]:
                    waveArr = []
                    print "   making waveArr: " + str(wavelength)
                    for obsNum in range(len(self.obsFileList)):
#                    for obsNum in [0,2]:
                        pixelArr = np.zeros((self.nxpix, self.nypix))
                        print "      obs: " + str(obsNum)
                        for x in range(self.nxpix):
                            for y in range(self.nypix):
                                pixelArr[x,y] = self.isGoodQE(obsNum, self.nxpix*y+x, wavelength)
                        waveArr.append(pixelArr)
                    dataQE.append(np.asarray(waveArr))

                print "DONE!"
                self.QEData = np.asarray(dataQE)
                np.save('QEData.npy',self.QEData)
                print " "

    def isGoodQE(self, obsNum, pixNum, wavelength,deadTime=100):
        #deadTime in microseconds for linearity correction. deadTime=0 applies no correction

        if (self.currentObs>=0):           
            try:
                self.obsFileList[self.currentObs]
                self.QEFileList[self.currentObs]
            except IndexError:
                print 'Missing observation or QE file'

            currentPix = self.selectedPixel
            currentObsNum = self.currentObs

            self.currentObs = obsNum
            self.selectedPixel = pixNum
            #self.goodPixelQELimit = 0.1

            QEFileObject = self.QEFileList[self.currentObs]
            ind = where(QEFileObject.data[:,0] == wavelength)[0][0]
            waveArr = QEFileObject.data[ind]   
#            print 'Obs: '+str(obsNum)+ ' pix: (' + str(self.selectedPixel/self.nxpix)+','+str(self.selectedPixel%self.nxpix)+')'+' wavelength: ' +str(wavelength)

            backGroundRate1, sb1 = self.getAverageTimeStream(self.currentObs,waveArr[-6],waveArr[-5])
#            print 'backGroundRate1: ' +str(backGroundRate1)
            backGroundRate2, sb2 = self.getAverageTimeStream(self.currentObs,waveArr[-2],waveArr[-1])
#            print 'backGroundRate2: ' +str(backGroundRate2)
            
            countRate, s = self.getAverageTimeStream(self.currentObs,waveArr[-4],waveArr[-3])
#            print 'countRate: ' +str(countRate)

            if deadTime>0:
                backGroundRate1=backGroundRate1/(1.0-backGroundRate1*deadTime*10.0**-6.0)
                backGroundRate2=backGroundRate2/(1.0-backGroundRate2*deadTime*10.0**-6.0)
                countRate=countRate/(1.0-countRate*deadTime*10.0**-6.0)
            
            backGroundRate = (backGroundRate1+backGroundRate2)/2.0
            countRate-=backGroundRate
            QE=QEFileObject.findQE(countRate,waveArr[0])

  
            #if background doesn't change significantly before and after
            #if (QE>0.01) & (abs(backGroundRate2 - backGroundRate1)<3*math.sqrt(backGroundRate)) & (countRate>3*math.sqrt(backGroundRate)):
            #if (abs(backGroundRate2 - backGroundRate1)<3*math.sqrt(backGroundRate)) & (countRate>1.0/2.0*math.sqrt(backGroundRate)):

            err=math.sqrt(backGroundRate)
            if backGroundRate<1:
                err=math.sqrt(5)
            if (abs(backGroundRate2 - backGroundRate1)<3*err):
#                print 'Obs: '+str(obsNum)+ ' pix: (' + str(self.selectedPixel/self.nxpix)+','+str(self.selectedPixel%self.nxpix)+')'+' wavelength: ' +str(wavelength)
#                print 'backGroundRate1: ' +str(backGroundRate1)
#                print 'backGroundRate2: ' +str(backGroundRate2)
#                print 'countRate: ' +str(countRate)
                self.selectedPixel = currentPix
                self.currentObs = currentObsNum 
                return QE
            else:
                self.selectedPixel = currentPix
                self.currentObs = currentObsNum 
                return -0.01


#            if QE > self.goodPixelQELimit:
#                return QE
#            else:
#                return -1


    def findGoodPixel(self):
        
        if (self.currentObs>=0):
            try:
                self.obsFileList[self.currentObs]
                self.QEFileList[self.currentObs]
            except IndexError:
                print 'Missing observation file'
            currentPix = self.selectedPixel
            goodPixels = []
            #self.goodPixelQELimit = 0.1

            QEFileObject = self.QEFileList[self.currentObs]
            wavelength = float(self.wv_combobox.currentText())
            print "Looking for good pixels in obs " +str(self.currentObs) +" at wavelength: " + str(wavelength)
            ind = where(QEFileObject.data[:,0] == wavelength)[0][0]
            waveArr = QEFileObject.data[ind]                                    #timing data for wavelength selected and current obs
            for x in range(self.nxpix):
                for y in range(self.nypix):
                    if self.isGoodQE(self.currentObs, self.nxpix*y+x, wavelength)>self.goodPixelQELimit:
                        goodPixels.append(self.nxpix*y+x)

#                    self.selectedPixel = self.nxpix*x+y
#                    backGroundRate1 = self.getAverageTimeStream(self.currentObs,waveArr[-6],waveArr[-5])
#                    backGroundRate2 = self.getAverageTimeStream(self.currentObs,waveArr[-2],waveArr[-1])
#                    backGroundRate = (backGroundRate1+backGroundRate2)/2.0
#                    countRate = self.getAverageTimeStream(self.currentObs,waveArr[-4],waveArr[-3])
#                    countRate-=backGroundRate
#                    QE=QEFileObject.findQE(countRate,waveArr[0])

#                    if QE > self.goodPixelQELimit:
#                        goodPixels.append(self.selectedPixel)
                    
            self.selectedPixel = currentPix
            self.good_pixel = goodPixels
            self.display_image()
            print "Found: " +str(self.pix_toTuple(self.good_pixel))

    def pix_toTuple(self, pixelArr):
        newArr = []
        for pix in pixelArr:
            newArr.append([pix%(self.nxpix),pix/self.nxpix])
        return newArr

    def start_pixel_select(self,event):
        #Mouse press returns x,y position of first pixel to be used in spectra
        startrawx,startrawy = event.pos().x(), event.pos().y()
        if hasattr(self,'scalex'):
            startpx = int(startrawx/self.scalex)
            startpy = int((self.nypix) - startrawy/self.scaley)
            startpix = self.nxpix*startpy+startpx
            print 'Pixel: ('+ str(startpx) +', '+ str(startpy)+')'
            self.pixel_label.setText('Pixel: ('+ str(startpx) +', '+ str(startpy)+')')
            self.selectedPixel=startpix
            self.display_image()
            self.plot_timestream()
            self.plot_QE()

    def getAverageTimeStream(self, obsNum, tstart, tend):
        x=self.selectedPixel%(self.nxpix)
        y=self.selectedPixel/self.nxpix


        #QEFileObject = self.QEFileList[obsNum]
        bmap = self.bmapList[obsNum]
        ti=int(tstart)
        tf=int(tend)
#        print '(t1,t2): ('+str(ti)+','+str(tf)+')'
        counts = zeros(tf-ti)
        h5file = openFile(str(self.obsFileList[obsNum]), 'r')

        for t in range(tf-ti):

            try:

                if bmap[y][x] == '':
                    continue
                try:
                    counts[t] += len(h5file.root._f_getChild(bmap[y][x])[ti+t])
                    if counts[t]<0:
                        counts[t]=0
                     
                except NoSuchNodeError:
                    counts[t]=0  

            except IndexError:
                print '(x,y): ' + '(', str(x)+','+str(y)+')'
                print 'pixel: ' + str(self.selectedPixel)
                print 'nxpix: ' + str(self.nxpix)
                print 'nypix: ' +str(self.nypix)
                print 'len bmap: ' + str(len(bmap))
                print 'len bmap[0]: ' + str(len(bmap[0])) 
                pass

        h5file.close()
#        if ti>1400:
#            print 'counts: ' +str(counts)

        #return np.mean(counts)

#        print 'std: ' + str(np.std(counts, ddof = 1))
        return np.median(counts), np.std(counts, ddof = 1)

    def plot_timestream(self):   
        print 'current obs: '+str(self.currentObs)
        obsfilename = self.obsFileList[self.currentObs]
        x=self.selectedPixel%(self.nxpix)
        y=self.selectedPixel/self.nxpix


#        i=self.startpy
#        j=self.startpx
        bmap = self.bmapList[self.currentObs]
#        self.ui.pixelpath.setText(str(bmap[i][j]))        
#        print "plotting time stream for" + str(bmap[i][j])
#        self.pixel_coord.setText(str(bmap[i][j]))
        print "plotting time stream for" + str(bmap[y][x])
        self.pixel_coord.setText(str(bmap[y][x]))
#        nphot=0

#        ti=self.startTime
        ti=0
        tf=self.obsTime
            
        if ti>tf:
            copyti = ti
            ti=tf
            tf=copyti
            print "WARNING: selected ti > tf. They were switched for you."
            
        counts = zeros(tf-ti)
        dcounts = zeros(tf-ti-1)
        h5file = openFile(str(obsfilename), 'r')
            


        for t in range(tf-ti):

            if bmap[y][x] == '':
                continue
            try:
                counts[t] += len(h5file.root._f_getChild(bmap[y][x])[ti+t])
                if t <= (tf-ti-2):
                    dcounts[t]= len(h5file.root._f_getChild(bmap[y][x])[ti+t+1])-len(h5file.root._f_getChild(bmap[y][x])[ti+t])
#                    if dcounts[t]<0:
#                        dcounts[t]=0
                if counts[t]<0:
                    counts[t]=0
                 
            except NoSuchNodeError:
                counts[t]=0  

        timesteps = xrange(ti,tf)
        dtimesteps = xrange(ti,tf-1)

        #print "plotting timestream of ", tf-ti, " seconds"
#        self.ui.countlabel.setText(str(sum(counts)))
        self.timestream_plot.canvas.ax.clear()
        self.timestream_plot.canvas.ax.plot(timesteps,counts)
        self.Dtimestream_plot.canvas.ax.clear()
        self.Dtimestream_plot.canvas.ax.plot(dtimesteps,dcounts)
        self.timestream_plot.format_labels('Time Stream', 'time (s)', 'counts')
        self.Dtimestream_plot.format_labels('D Time Stream', 'time (s)', 'dcounts/dt')

        try:
            self.plot_timeAreas()
        except IndexError:
            pass
    

        self.timestream_plot.canvas.draw()
        self.Dtimestream_plot.canvas.draw()

        h5file.close()
        print "done plotting timestream"

    def addNewObsFile(self):
        newdatafile = QFileDialog.getOpenFileName(parent=None, caption=QString(str("Choose Obs File")),directory = "/Scratch/QEData/",filter=QString(str("H5 (*.h5)")))
        self.initializeObsFile(newdatafile)

    def initializeObsFile(self, fn):
        newdatafile = fn
        if len(newdatafile)!=0:
            self.obsFileList.append(newdatafile)
            self.currentObs = self.numObservations
            self.beammapfile = str(newdatafile)
            #self.beammapfile = '/Scratch/LABTEST/20140910/beammap_SCI6_B140731-Force_20140910.h5'
            self.loadbeammap()
            self.bmapList.append(self.bmap)
#            self.currentObs = len(self.obsFileList)-1
#            try:
#                self.QEFileList[self.currentObs]
#            except IndexError:
#                i=len(self.QEFileList)
#                while i <= self.currentObs:
#                    self.QEFileList.append(None)
#                    i+=1
            print 'current obs: '+str(self.currentObs)
            self.obsNum_label.setText('Measurement: ' + str(self.currentObs))
            splitname = newdatafile.split("/")
            justname = splitname[-1]
            
            self.obsfile_label.setText(str(justname))

            self.blankscene = QGraphicsScene()
            self.blankscene.clear()
            self.tv_image.setScene(self.blankscene)
            self.tv_image.show()
#            self.ui.histogram_plot.canvas.ax.clear()
#            self.ui.histogram_plot.canvas.format_labels()
#            self.ui.histogram_plot.canvas.draw()

            self.get_ut()
            
            print 'Making image'
            fn = self.obsFileList[self.currentObs].split("/")[-1].split(".")[0]
            self.imagefile = "./TV_frame_"+str(fn)+".png"
            if isfile(self.imagefile):
                self.display_image()
            else:
                self.make_image()


    def addNewQEFile(self):
        newQEfile=QFileDialog.getOpenFileName(parent=None, caption=QString(str("Choose Obs File")),directory = "/Scratch/QEData/",filter=QString(str("txt (*.txt)")))
        angle = self.askForAngle()
        if angle != None:
            self.initializeQEFile(newQEfile,angle)

    def askForAngle(self):
        angle, ok = QInputDialog.getInt(self,'Input Dialog here', 'Enter angle to nearest degree:')
        if ok:
            return angle
        else:
            print 'Cancelled'
            return None

    def initializeQEFile(self, fn, angle):
        newQEfile = fn
        if len(newQEfile) != 0:
#            angle, ok = QInputDialog.getInt(self,'Input Dialog here', 'Enter angle to nearest degree:')
#            if ok:
#                self.QEFileList.append(QEfile(newQEfile,angle))
#                self.currentObs = self.numObservations
#                self.currentObs = len(self.QEFileList)-1
#                print len(self.QEFileList)
#                print 'current obs: '+str(self.currentObs)
#                self.obsNum_label.setText('Measurement: ' + str(self.currentObs))
#                splitname = newQEfile.split("/")
#                justname = splitname[-1]
#                self.QEfile_label.setText(str(justname))
#               self.addQETiming()
#                self.populateTimeAreaData()
#            else:
#                print 'cancelled'
                

            self.QEFileList.append(QEfile(newQEfile,angle))
            self.currentObs = self.numObservations
#                self.currentObs = len(self.QEFileList)-1
#            print len(self.QEFileList)
            print 'current obs: '+str(self.currentObs)
            self.obsNum_label.setText('Measurement: ' + str(self.currentObs))
            splitname = newQEfile.split("/")
            justname = splitname[-1]
            self.QEfile_label.setText(str(justname))
#               self.addQETiming()
            self.populateTimeAreaData()

#            try:
#                self.obsFileList[self.currentObs]
#            except IndexError:
#                i=len(self.obsFileList)
#                while i <= self.currentObs:
#                    self.obsFileList.append(None)
#                    i+=1


    def addQETiming(self):
        if (self.currentObs >= 0) & (len(self.QEFileList) >= (self.currentObs+1)):
            
            QEdata=self.QEFileList[self.currentObs].data
            
            for ind in range(len(QEdata)):
                self.addWaveToFrame(QEdata[ind],ind+1)

    def make_image(self):
        obsfilename = self.obsFileList[self.currentObs]
        bmap = self.bmapList[self.currentObs]
        nphot=0
#        ti=self.startTime
#        tf=self.endTime
#        if ti>tf:
#            copyti = ti
#            ti=tf
#            tf=copyti
#            print "WARNING: selected ti > tf. They were switched for you."

        h5file = openFile(str(obsfilename), 'r')



        
        all_photons = []
        for j in range(self.nxpix*self.nypix):
            all_photons.append([])
        
        counts = zeros((self.nypix,self.nxpix))
        
        for i in xrange(self.nypix):
            for j in xrange(self.nxpix):
                    if bmap[i][j] == '':
                        counts[i][j]=0
                        continue
                    try:
#                        print '('+str(i)+', '+str(j)+')'
                        photons = concatenate(h5file.root._f_getChild(bmap[i][j])[0:self.obsTime])
#                        photons = concatenate(h5file.root._f_getChild(bmap[i][j])[0:3])
                        counts[i][j]=len(photons)

                        if counts[i][j]<0:
                            counts[i][j]=0

                        
                    except NoSuchNodeError:
                        counts[i][j]=0

                    nphot += counts[i][j]
    
        photon_count = counts
        im = photon_count
        
        photon_count = flipud(photon_count)
#        if self.ui.man_contrast.isChecked():
#            self.vmax = self.ui.vmax.value()
#            self.vmin = self.ui.vmin.value()
#        else:
#            indices = sort(reshape(photon_count,self.nxpix*self.nypix))
#            self.satpercent=2.5
#            brightest = int((self.satpercent/100.0)*(self.nxpix*self.nypix))
#            self.vmax = indices[(-1*brightest)]
#            self.vmin = 0

        indices = sort(reshape(photon_count,self.nxpix*self.nypix))
        self.satpercent=2.5
        brightest = int((self.satpercent/100.0)*(self.nxpix*self.nypix))
        self.vmax = indices[(-1*brightest)]
        self.vmin = 0
            
        fig = plt.figure(figsize=(0.01*self.nxpix,0.01*self.nypix), dpi=100, frameon=False)
        im = plt.figimage(photon_count, cmap='gray', vmin = self.vmin, vmax = self.vmax)
        fn = self.obsFileList[self.currentObs].split("/")[-1].split(".")[0]
        self.imfile = "TV_frame_"+str(fn)+".png"
        plt.savefig(self.imfile, pad_inches=0)
        print "done making image."
            
        self.display_image()
        h5file.close()

    def display_image(self):
        #search directory for image
        fn = self.obsFileList[self.currentObs].split("/")[-1].split(".")[0]
        self.imagefile = "./TV_frame_"+str(fn)+".png"
        
        self.tv_image.setGeometry(self.position.x(), self.position.y()-8, self.scalex*(self.nxpix)+4, self.scaley*(self.nypix)+4)
        if isfile(self.imagefile):
            pix = self.makepixmap(self.imagefile, scalex=self.scalex, scaley=self.scaley)
            #display pixmap
            self.scene = QGraphicsScene()
            self.scene.addPixmap(pix)
            
#            if self.histogram_pixel != []:
#                for p in self.histogram_pixel:
#                    x = p%(self.nxpix)
#                    y = (self.nypix-1) - p/self.nxpix                    
#                    self.scene.addRect(self.scalex*(x),self.scaley*(y),(self.scalex),(self.scaley), Qt.blue)
            if self.selectedPixel != None:
                x=self.selectedPixel%(self.nxpix)
                y = (self.nypix-1) - self.selectedPixel/self.nxpix
                self.scene.addRect(self.scalex*(x),self.scaley*(y),(self.scalex),(self.scaley), Qt.blue)
            for pix in self.good_pixel:
                x=pix%(self.nxpix)
                y = (self.nypix-1) - pix/self.nxpix
                self.scene.addRect(self.scalex*(x),self.scaley*(y),(self.scalex),(self.scaley), Qt.green)
                if pix == self.selectedPixel:
                    self.scene.addRect(self.scalex*(x),self.scaley*(y),(self.scalex),(self.scaley), Qt.darkCyan)
                    
            self.tv_image.setScene(self.scene)
            self.tv_image.show()
            #os.remove(str(imagefile)) #want to keep this around for saving purposes
        else: 
            print 'No image found: ' +str(self.imagefile)
            self.blankscene = QGraphicsScene()
            self.blankscene.clear()
            self.tv_image.setScene(self.blankscene)
            self.tv_image.show()

    def get_ut(self):
        obsfilename = self.obsFileList[self.currentObs]
        h5file = openFile(str(obsfilename), mode = "r")
        htable = h5file.root.header.header.read()
        self.obsTime = int(htable["exptime"])
        h5address = h5file.root.beammap.beamimage.read()[0][0]
        h5time = int(h5address.split('t')[1])
        try:
            self.ut = int(htable["unixtime"])
        except ValueError:
            print "unixtime not found, checking for deprecated ut field" 
            self.ut = int(htable["ut"])

        if self.ut != h5time:
            self.ut = h5time
            
        h5file.close()
        self.bmap = self.bmapList[self.currentObs]
        for i in xrange(self.nypix):
            for j in xrange(self.nxpix):
                head = str(self.bmap[i][j])
                if head.split('t')[0] == '':
                    self.bmap[i][j] = head + 't' + str(self.ut)
                else:
                    self.bmap[i][j] = head.split('t')[0] + 't' + str(self.ut)
        

        print "Pixel addresses updated in beammap"

    def makepixmap(self, imagefile, scalex=1, scaley=1):
        '''Given an image file this function converts them to pixmaps to be displayed by QT gui'''
        qtimage = QImage(imagefile)
        width, height = qtimage.size().width(), qtimage.size().height()
        qtimage = qtimage.scaled(width*scalex,height*scaley)
        pix = QPixmap.fromImage(qtimage)
        return pix

    def loadbeammap(self):
        bmfile = openFile(self.beammapfile, 'r')
        #read beammap in to memory to create beam image
        self.bmap = bmfile.root.beammap.beamimage.read()
        self.nxpix = shape(self.bmap)[1]
        self.nypix = shape(self.bmap)[0]
        bmfile.close()
        print "Beammap loaded from " +str(self.beammapfile)

        self.scalex=float(self.imagex/float(self.nxpix))
        self.scaley=float(self.imagey/float(self.nypix))

    def addtomainframe(self):
        self.button_0 = QPushButton("button 0")
        self.vbox.addWidget(self.button_0)
        self.main_frame.setLayout(self.vbox)
        self.setCentralWidget(self.main_frame)

    def plot_timeAreas(self):
        QEFileObject = self.QEFileList[self.currentObs]
        for waveArr in QEFileObject.data:
#            print waveArr
            self.timestream_plot.canvas.ax.axvspan(waveArr[5], waveArr[6], facecolor='r', alpha=0.5)
            self.timestream_plot.canvas.ax.axvspan(waveArr[7], waveArr[8], facecolor='g', alpha=0.5)
            self.timestream_plot.canvas.ax.axvspan(waveArr[9], waveArr[10], facecolor='r', alpha=0.5)
            self.Dtimestream_plot.canvas.ax.axvspan(waveArr[5], waveArr[6], facecolor='r', alpha=0.5)
            self.Dtimestream_plot.canvas.ax.axvspan(waveArr[7], waveArr[8], facecolor='g', alpha=0.5)
            self.Dtimestream_plot.canvas.ax.axvspan(waveArr[9], waveArr[10], facecolor='r', alpha=0.5)
        self.timestream_plot.format_labels('Time Stream', 'time (s)', 'counts')
        self.Dtimestream_plot.format_labels('D Time Stream', 'time (s)', 'dcounts/dt')
        self.timestream_plot.canvas.draw()
        self.Dtimestream_plot.canvas.draw()


    def timeChanged(self):
        source = self.sender()
        value = int(source.value())
        self.QEFileList[self.currentObs].data[int(source.y), int(source.x)]=value
        self.plot_timestream()
        self.plot_QE()

    def addWaveToFrame(self):
        #QE calibration
        #self.test_button = QPushButton('test')
        #QObject.connect(self.test_button, SIGNAL('clicked()'),self.toggleEnable)
        obsfile_button = QPushButton("upload obs file")
        QObject.connect(obsfile_button, SIGNAL('clicked()'),self.addNewObsFile)
        #if obs file exists....
        self.obsfile_label = QLabel("")
        self.obsfile_label.setFrameShape(QFrame.Box)
        self.obsfile_label.setMinimumWidth(300)
        self.obsfile_label.setMaximumHeight(25)
        QEfile_button = QPushButton("upload QE file")
        QObject.connect(QEfile_button, SIGNAL('clicked()'),self.addNewQEFile)
        #if obs file exists....
        self.QEfile_label = QLabel("")
        self.QEfile_label.setMinimumWidth(300)
        self.QEfile_label.setMaximumHeight(25)
        self.QEfile_label.setFrameShape(QFrame.Box)

        



        #if obs file exists.... populateTimeAreaData()

        self.vbox_cal = QVBoxLayout()
        hbox_obs = QHBoxLayout()
        hbox_obs.addWidget(obsfile_button)
        hbox_obs.addWidget(self.obsfile_label)
        self.vbox_cal.addLayout(hbox_obs)
        hbox_QE = QHBoxLayout()
        hbox_QE.addWidget(QEfile_button)
        hbox_QE.addWidget(self.QEfile_label)
        self.vbox_cal.addLayout(hbox_QE)

        self.cal_QFrame.setLayout(self.vbox_cal)


    def populateTimeAreaData(self):

        wave_label = QLabel("wavelength")
        t1_label = QLabel("Trough 1 Start")
        t2_label = QLabel("Trough 1 End")
        f1_label = QLabel("Flux Start")
        f2_label = QLabel("Flux End")
        t3_label = QLabel("Trough 2 Start")
        t4_label = QLabel("Trough 2 End")

        grid_wave = QGridLayout()
        grid_wave.addWidget(wave_label,0,0)
        grid_wave.addWidget(t1_label,0,1)
        grid_wave.addWidget(t2_label,0,2)
        grid_wave.addWidget(f1_label,0,3)
        grid_wave.addWidget(f2_label,0,4)
        grid_wave.addWidget(t3_label,0,5)
        grid_wave.addWidget(t4_label,0,6)


        QEFileObject = self.QEFileList[self.currentObs]
        row=0
        for waveArr in QEFileObject.data:
            row+=1
            wave_label = QLabel(str(waveArr[0]))
#            t1_text = QLineEdit(str(waveArr[5]))
            t1_spinbox = QSpinBox()
            t1_spinbox.setMaximum(100000)
            t1_spinbox.setValue(int(waveArr[5]))
            t1_spinbox.x = 5
            t1_spinbox.y = row-1
            QObject.connect(t1_spinbox, SIGNAL('editingFinished()'),self.timeChanged)
            t2_spinbox = QSpinBox()
            t2_spinbox.setMaximum(100000)
            t2_spinbox.setValue(int(waveArr[6]))
            t2_spinbox.x = 6
            t2_spinbox.y = row-1
            QObject.connect(t2_spinbox, SIGNAL('editingFinished()'),self.timeChanged)
            f1_spinbox = QSpinBox()
            f1_spinbox.setMaximum(100000)
            f1_spinbox.setValue(int(waveArr[7]))
            f1_spinbox.x = 7
            f1_spinbox.y = row-1
            QObject.connect(f1_spinbox, SIGNAL('editingFinished()'),self.timeChanged)
            f2_spinbox = QSpinBox()
            f2_spinbox.setMaximum(100000)
            f2_spinbox.setValue(int(waveArr[8]))
            f2_spinbox.x = 8
            f2_spinbox.y = row-1
            QObject.connect(f2_spinbox, SIGNAL('editingFinished()'),self.timeChanged)
            t3_spinbox = QSpinBox()
            t3_spinbox.setMaximum(100000)
            t3_spinbox.setValue(int(waveArr[9]))
            t3_spinbox.x = 9
            t3_spinbox.y = row-1
            QObject.connect(t3_spinbox, SIGNAL('editingFinished()'),self.timeChanged)
            t4_spinbox = QSpinBox()
            t4_spinbox.setMaximum(100000)
            t4_spinbox.setValue(int(waveArr[10]))
            t4_spinbox.x = 10
            t4_spinbox.y = row-1
            QObject.connect(t4_spinbox, SIGNAL('editingFinished()'),self.timeChanged)

            grid_wave.addWidget(wave_label,row,0)
            grid_wave.addWidget(t1_spinbox,row,1)
            grid_wave.addWidget(t2_spinbox,row,2)
            grid_wave.addWidget(f1_spinbox,row,3)
            grid_wave.addWidget(f2_spinbox,row,4)
            grid_wave.addWidget(t3_spinbox,row,5)
            grid_wave.addWidget(t4_spinbox,row,6)


        angle = QEFileObject.angle
        self.angle_label = QLabel('Angle: '+str(angle))
        
        doneCal_button = QPushButton('Done Calibration')
        QObject.connect(doneCal_button, SIGNAL('clicked()'),self.doneCal)
        self.vbox_cal.addWidget(self.angle_label)
        self.vbox_cal.addWidget(doneCal_button)

        self.vbox_cal.addLayout(grid_wave)

    def switchObs(self):
        source = self.sender()
        self.currentObs = int(source.num)
        self.obsNum_label.setText('Measurement: ' + str(self.currentObs))
        self.display_image()
        if self.selectedPixel != None:
            self.plot_timestream()

    def plot_QE(self):
        if (len(self.QEFileList) > 0) & (len(self.obsFileList) == len(self.QEFileList)) & (self.selectedPixel != None):
            self.wvQE_plot.canvas.ax.clear()
            self.angleQE_plot.canvas.ax.clear()
            wvQEArr = []
            angleQEArr = []
            for obsNum in range(len(self.obsFileList)):
                
                QEFileObject = self.QEFileList[obsNum]
                wvQEArr.append(zeros((len(QEFileObject.data),2)))
                
                for wvNum in range(len(QEFileObject.data)):
                    waveArr = QEFileObject.data[wvNum]
                    backGroundRate1, sb1 = self.getAverageTimeStream(obsNum,waveArr[-6],waveArr[-5])
                    backGroundRate2, sb2 = self.getAverageTimeStream(obsNum,waveArr[-2],waveArr[-1])
                    backGroundRate = (backGroundRate1+backGroundRate2)/2.0
                    countRate, s = self.getAverageTimeStream(obsNum,waveArr[-4],waveArr[-3])
                    countRate-=backGroundRate
                    #print "obsnum: " + str(obsNum)
                    #print "wavelength: " +str(waveArr[0])
                    #print "backGroundRate: " +str(backGroundRate)
                    #print "countRate: " +str(countRate)
                    
                    QE=QEFileObject.findQE(countRate,waveArr[0])
                    #print "QE: "+str(QE)
                    if QE < 0:
                        QE = 0
#                    print 'QE: ' + str(QE)+' wavelength: '+str(waveArr[0])
                    wvQEArr[obsNum][wvNum,:]=[waveArr[0], QE]

                    if int(waveArr[0]) == int(self.wv_combobox.currentText()):
                        angleQEArr.append([QEFileObject.angle, QE])

                self.wvQE_plot.canvas.ax.plot(wvQEArr[obsNum][:,0],wvQEArr[obsNum][:,1], label = str(obsNum))             

            angleQEArr=np.reshape(angleQEArr,(len(angleQEArr),2))
#            self.angleQE_plot.canvas.ax.plot(angleQEArr[:,0], angleQEArr[:,1])
            self.angleQE_plot.canvas.ax.plot(angleQEArr[:,0], angleQEArr[:,1],'bo')
#            print 'AngleArr: ' + str(angleQEArr)
#            print 'x: ' +str(angleQEArr[:,0])
#            print 'y: ' +str(angleQEArr[:,1])

#            self.wvQE_plot.canvas.ax.legend(loc='upper right')
#            handles, labels = self.wvQE_plot.canvas.ax.get_legend_handles_labels()
#            self.wvQE_plot.canvas.ax.legend(handles, labels, loc='upper right')
            self.wvQE_plot.canvas.ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
            self.angleQE_plot.format_labels('Polarization', 'Angle (degrees)', 'QE')
            self.wvQE_plot.format_labels('QE', 'wavelength (nm)', 'QE')
            self.angleQE_plot.canvas.draw()
            self.wvQE_plot.canvas.draw()

        

    def doneCal(self):

        QEFileExists = False
        obsFileExists = False
        try:
            QEFileExists = isfile(str(self.QEFileList[self.currentObs].filename))
        except IndexError:
            pass
        try:
            obsFileExists = isfile(str(self.obsFileList[self.currentObs]))
        except IndexError:
            pass
        
        if QEFileExists & (not obsFileExists):
            self.QEFileList.pop(self.currentObs)
#            self.currentObs-=1
            self.currentObs = self.numObservations-1
        elif obsFileExists & (not QEFileExists):
            self.obsFileList.pop(self.currentObs)
#            self.currentObs-=1
            self.currentObs = self.numObservations-1
        elif obsFileExists & QEFileExists:
            self.numObservations+=1
            self.QEFileList[self.currentObs].saveParam()
            obsNum_button = QPushButton('('+str(self.currentObs)+')')
            obsNum_button.num = self.currentObs
            QObject.connect(obsNum_button, SIGNAL('clicked()'),self.switchObs)
            editObs_button = QPushButton('Edit')
            editObs_button.num = self.currentObs
            angle_label_temp = QLabel(str(self.angle_label.text()))
            obsfile_label_temp = QLabel(str(self.obsfile_label.text()))
            QEfile_label_temp = QLabel(str(self.QEfile_label.text()))
            hbox0 = QHBoxLayout()
            vbox0 = QVBoxLayout()
            vbox0.addWidget(angle_label_temp)
            vbox0.addWidget(obsfile_label_temp)
            vbox0.addWidget(QEfile_label_temp)
            hbox0.addWidget(obsNum_button)
            hbox0.addLayout(vbox0)
            hbox0.addWidget(editObs_button)
            self.vbox_file.insertLayout(self.currentObs, hbox0)

        self.cal_QFrame.hide()
        self.vbox_QEcal.removeWidget(self.cal_QFrame)
        self.cal_QFrame = QFrame()
        self.cal_QFrame.hide()
        self.vbox_QEcal.addWidget(self.cal_QFrame)
        self.file_QFrame.show()
        

    def toggleEnable(self):
#        self.grid_wave.hide()
        for i in range(7):
            if self.grid_wave.itemAtPosition(1,i).widget().isVisible():
                self.grid_wave.itemAtPosition(1,i).widget().hide()
            else:
                self.grid_wave.itemAtPosition(1,i).widget().show()

    def addMeasurementFile(self):
        self.file_QFrame.hide()
        self.addWaveToFrame()
        self.cal_QFrame.show()


    def create_main_frame(self):
        self.main_frame = QWidget()
        self.tv_image = QGraphicsView()
        self.tv_image.setMinimumSize(self.imagex,self.imagey)
        self.tv_image.setMaximumSize(self.imagex,self.imagey)
        self.tv_image.setMouseTracking(True)
        self.tv_image.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.tv_image.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.tv_image.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignTop)
        self.pixel_label = QLabel("pixel")
        self.obsNum_label = QLabel('Measurement: ' + str(self.currentObs))
        self.pixel_coord = QLabel('pixel coord')
        self.pixel_button = QPushButton('Toggle good pixel')
        QObject.connect(self.pixel_button, SIGNAL('clicked()'),self.toggleGoodPixel)
        self.goodpixel_button = QPushButton('Find good pixels')
        QObject.connect(self.goodpixel_button, SIGNAL('clicked()'),self.findGoodPixel)
        self.timestream_plot = MPL_Widget()
        self.timestream_plot.format_labels('Time Stream', 'time (s)', 'counts')
        self.Dtimestream_plot = MPL_Widget()
        self.Dtimestream_plot.format_labels('D Time Stream', 'time (s)', 'dcounts/dt')

        self.file_QFrame = QFrame()
        self.cal_QFrame = QFrame()                       #calibration box
        self.cal_QFrame.hide()
        addMeasurement_button = QPushButton('Add a new Measurement')
        QObject.connect(addMeasurement_button, SIGNAL('clicked()'),self.addMeasurementFile)


        self.wvQE_plot = MPL_Widget()
        self.wvQE_plot.format_labels('QE', 'wavelength (nm)', 'QE')
        self.wv_combobox = QComboBox()
        self.wv_combobox.addItems(['400', '450', '500', '550', '600', '650', '700', '750', '800', '850', '900', '950', '1000', '1050', '1100', '1150', '1200', '1250', '1300'])
        QObject.connect(self.wv_combobox, SIGNAL('currentIndexChanged(QString)'),self.plot_QE)
        self.average_button = QPushButton('Find Obs Average')
        QObject.connect(self.average_button, SIGNAL('clicked()'),self.getAvgQE)
        self.compare_button = QPushButton('Compare Averages')
        QObject.connect(self.compare_button, SIGNAL('clicked()'),self.compareQE)
        self.histogram_button = QPushButton('Make Histograms')
        QObject.connect(self.histogram_button, SIGNAL('clicked()'),self.histogramQE)
        self.angleQE_plot = MPL_Widget()
        self.angleQE_plot.format_labels('Polarization', 'Angle (degrees)', 'QE')
        



        hbox = QHBoxLayout()                        #main box

        vbox_plot = QVBoxLayout()                        #plot box
        hbox_tv = QHBoxLayout()
        hbox_tv.addWidget(self.tv_image)
        vbox_pixellabels = QVBoxLayout()
        vbox_pixellabels.addWidget(self.pixel_label)
        vbox_pixellabels.addWidget(self.obsNum_label)
        vbox_pixellabels.addWidget(self.pixel_coord)
        vbox_pixellabels.addWidget(self.pixel_button)
        vbox_pixellabels.addWidget(self.goodpixel_button)
        hbox_tv.addLayout(vbox_pixellabels)
        vbox_plot.addLayout(hbox_tv)
        vbox_plot.addWidget(self.timestream_plot)
        vbox_plot.addWidget(self.Dtimestream_plot)


        self.vbox_file = QVBoxLayout()                   #add New Files box
        self.vbox_file.addWidget(addMeasurement_button)
        self.file_QFrame.setLayout(self.vbox_file)

        self.vbox_QEcal = QVBoxLayout()
        self.vbox_QEcal.addWidget(self.file_QFrame)
        self.vbox_QEcal.addWidget(self.cal_QFrame)


        vbox_QEplot = QVBoxLayout()
        vbox_QEplot.addWidget(self.wvQE_plot)
        hbox_QEbuttons = QHBoxLayout()
        hbox_QEbuttons.addWidget(self.wv_combobox)
        hbox_QEbuttons.addWidget(self.average_button)
        hbox_QEbuttons.addWidget(self.compare_button)
        hbox_QEbuttons.addWidget(self.histogram_button)
        vbox_QEplot.addLayout(hbox_QEbuttons)
        vbox_QEplot.addWidget(self.angleQE_plot)


        hbox.addLayout(vbox_plot)
        hbox.addLayout(self.vbox_QEcal)
        hbox.addLayout(vbox_QEplot)


        self.main_frame.setLayout(hbox)
        self.setCentralWidget(self.main_frame)

    def create_status_bar(self):
        self.status_text = QLabel("Awaiting orders.")
        self.statusBar().addWidget(self.status_text, 1)

    def on_about(self):
        msg = """ Message to user goes here.
        """
        QMessageBox.about(self, "MKID-ROACH software demo", msg.strip())

    def closeEvent(self, event):
#        for obsNum in range(len(self.obsFileList)):
#            self.imagefile = "./TV_frame"+str(self.currentObs)+".png"
#            if isfile(self.imagefile):
#                os.remove(str(self.imagefile))
#                print 'Removed: ' +str(self.imagefile)
        time.sleep(.1)
        QtCore.QCoreApplication.instance().quit

if __name__ == "__main__":
    app = QApplication(sys.argv)
    myapp = QECalibration()
    myapp.show()
    app.exec_()


