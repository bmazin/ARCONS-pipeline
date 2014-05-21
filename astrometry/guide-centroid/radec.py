'''
To Do: 
1.write down function for converting guide pixel to arcons pixel and use the same function for calculating its inverse to convert arcons back to guide pixel
2.Implement a better way to find reference guide pixel and reference arcons pixel
'''

from functions import *
import pyfits
import numpy as np
from astropy import wcs
from scipy.optimize import fsolve

inputFile = 'runRecord.txt'  

class radec(object):

    def __init__(self,recordFile='runRecord.txt',tolError=400,caldir='./cal/',refPixGuide=[548,370],refPixARCONS=[11.8688,12.9538]):
    
        self.sortList = []
        #notice that self.timeStamp contains timeStamp in seconds(integers)
        self.timeStampList = []
        self.caldir = caldir
        
        #initialize ARCONS parameters      
        self.guidePlate=3.5833*10**(-5)
        self.arconsPlate=1.208*10**(-4)
        self.dimX=44
        self.dimY=46
        self.refPixGuide = refPixGuide
        self.refPixARCONS = refPixARCONS
        self.sufix = ''
        
        #sorting the list of files in time-ascending order and discard any images that have error greater than 400 pixels
        fileList = []
        recordFile = open(recordFile,'r')
        while True:
            line = recordFile.readline()
        
            if line == '':
                break
            fileName,error = line.split()
            error = float(error)    
            
            
            if error < tolError:
                fileList.append(fileName)

        tempList1 = []
        tempList2 = []
        self.suffix = fileList[6:]
        
        for nfile in fileList:
            timeStamp = nfile[0:6]
            sec = timeConvert(timeStamp)
            tempList1.append(sec)
        
        #time-ascending order
        tempList1.sort() 

        for ntime in tempList1:
            timeStr = timeConvert(ntime)
            outputFileName = timeStr + suffix
            tempList2.append(outputFileName)        
        
        #output sorted file and timeStamp lists
        self.sortList = tempList2
        self.timeStampList = tempList1

    def centroid(self,worldCoor,refPixGuide=[548,370],refPixARCONS=[11.8688,12.9538]):
        '''
        Convert world coordinate into arcons coordinate.
        '''
        guidePixel = [] 

        for nfile in self.sortList:
            #include the prefix for file directory
            nfiledir = self.caldir + nfile
            try:
                pix = _world2pix(nfiledir,[worldCoor],'params.npz')
            except:
                pix = _world2pix(nfiledir,[worldCoor], None)
            guidePixel.append([nfile[0:6],pix[0]]) 

       
        x = np.array([coorX for coorX in (guidePixel[index][1][0] for index in range(len(guidePixel)))])
        y = np.array([coorY for coorY in (guidePixel[index][1][1] for index in range(len(guidePixel)))])
        
        #difference in X,Y ARCONS coordinate relative to the reference pixel
        deltaX = (x-self.refPixGuide[0])*self.guidePlate/self.arconsPlate
        deltaY = (y-self.refPixGuide[1])*self.guidePlate/self.arconsPlate
        
        #the reason for the minux sign in deltaY is because the axis is inverted
        arconsX = self.refPixARCONS[0] + deltaX
        arconsY = self.refPixARCONS[1] - deltaY

        index = 0
        arconsCoor = []
        print arconsX
        #eliminate any coordinate that is out of the range of the plate (0<x<44,0<y<46)
        for coor1 in arconsX:
            if coor1 > 0 and coor1 < 44 and arconsY[index] > 0 and arconsY[index] < 46:
                arconsCoor.append([coor1,arconsY[index],self.sortList[index][0:6]])
            index += 1
        return arconsCoor

    def photonMapping(self,timeStamp,xPhotonPixel,yPhotonPixel):
        #the timestamp variable has to be a string input of format eg: 041527
        
        
        #first convert timeStamp to seconds in interger and compare with the list to find the closest neighbor
        timeStamp = timeConvert(timeStamp)
        
        #convert to np array for faster processing
        timeStampListNp = np.array[self.timeStampList]
        
        def _find_nearest(array,value):
        #find closest neighbor from a given array for a given value
            idx = (np.abs(array-value)).argmin()
            return array[idx]
        
        matchedTime = _find_nearest(timeStampListNp,timeStamp)
        #file name of the closest matching time stamp
        matchedFile = timeConvert(matchedTime) + self.suffix
        
        deltaX = xPhotonPixel - self.refPixARCONS[0]
        deltaY = yPhotonPixel + self.refPixARCONS[1]
        
        #find (x,y) in guide pixel coordinate
        x = deltaX*(self.arconsPlate/self.guidePlate) + self.refPixGuide[0]
        y = deltaY*(self.arconsPlate/self.guidePlate) + self.refPixGuide[1]
        
        filePath = self.caldir + matchedFile    
        
        #return RA,DEC world coordinate
        world = _pix2world(filePath,[[x,y]])
        
        #RA,DEC in degreess
        RA = world[0][0]
        DEC = world[0][1]
        
        return RA,DEC

if __name__ == '__main__':
 
    test = radec(tolError=10000)
    
    nlist = test.centroid(worldCoor=[98.172422,-0.031578365])
    print nlist

'''
coor = [[14.02,-0.18],[14.02,-0.18],[14.02,-0.18],[14.02,-0.18],[14.02,-0.18],[14.02,-0.18],[8.08,11.57],[8.52,11.7],[8.9,11.9],[8.96,11.467],[22.1,15.67],[21.7,15.55],[21.39,15.05],[21.01,15.42],[21.88,14.93],[35.36,18.8],[36,18.88],[35.85,19.38],[3.7,24.45],[3.83,24.2],[3.58,24.45],[4.07,24.577],[4.07,24.7],[16.9,27.91],[17.3,27.42],[16.69,27.3],[16.8,27.05],[16.07,27.916],[29.55,31],[29.18,31.62],[29.3,31.74],[29.18,31.74],[29.3,31.37],[9.14,39.787],[9.27,39.416],[9.148,39.41],[8.9,39.1688],[8.28,40.03],[8.28,40.03],[22,43.86],[21.63,43.62],[21.88,43.25],[21.637,43.37],[21.143,44.3624]]

d = lambda c1,c2: np.sqrt((c1[0]-c2[0])**2+(c1[1]-c2[1])**2) 

index=0
dif = []
for item1 in nlist:
	dif.append(d(item1,coor[index]))
	index += 1

avg = sum(dif)/len(dif)

print 'average = %s' %avg


dev = []
for item in dif:
	dev.append((item-avg)**2)

stddev = np.sqrt(sum(dev)/len(dev))

print 'standar deviation = %s' %stddev
'''




		
