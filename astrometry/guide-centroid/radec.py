'''
To Do
1.write down function for converting guide pixel to arcons pixel and use the same function for calculating its inverse to convert arcons back to guide pixel
2.Implement a better way to find reference guide pixel and reference arcons pixel
'''

from functions import *
import pyfits
import numpy as np
from astropy import wcs
from scipy.optimize import fsolve
import itertools
import h5py
  

class radec(object):

    def __init__(self,tolErr=2500,caldir='/Scratch/larry/psr0614/',refPixGuide=[565,366],refPixARCONS=[11.18,13.5],suffix='_offCal_rotCal.fits'):
    
        self.sortList = []
        #notice that self.timeStamp contains timeStamp in seconds(integers)
        self.timeStampList = []
        self.caldir = caldir
        self.tolErr = tolErr
        #record the sorted file and corresponding directory
        self.dictionary = {}
        
        #initialize ARCONS parameters      
        self.guidePlate = 3.5833*10**(-5)
        self.arconsPlate = 1.208*10**(-4)
        self.dimX = 44
        self.dimY = 46
        self.refPixGuide = refPixGuide
        self.refPixARCONS = refPixARCONS
        self.suffix = suffix
        
        #sorting the list of files in time-ascending order and discard any images that have error greater than 2500 pixel distance in default
        
        #first, find all the directories with different date in the main calibration directory
        directories = []
        for subdir in os.listdir(self.caldir):
            directories.append(self.caldir+subdir+'/')
        
        #then sort the files in each directory separately. Create a dictionary and have the directory name to be keyword and sorted list be the argument
        for directory in directories:
            calList = []
            for calFile in os.listdir(directory):
                #firstly, fix this annoying error. Namely, the ambiguity in nameing the file (063960=064000!!) which causes error
                if calFile[4:6] == '60':
                    timeStamp = calFile[0:6]
                    #this double conversion will convert to the right format (i.e no 60 present)
                    newFormat = timeConvert(timeConvert(timeStamp)) + self.suffix
                    newFileName = directory + newFormat
                    oldFileName = directory + calFile
                    os.rename(oldFileName,newFileName)
                    calFile = newFormat
                    
                #skip loop if it is not a fits file
                if calFile[-5:] != '.fits':
                    continue
                imageList = pyfits.open(directory+calFile)
                header = imageList[0].header
                if header['CALERR'] < self.tolErr:
                    calList.append(calFile[0:6])
                imageList.close()
            #append dictionary entry. Here, the calList contains only the time stamp of each file(in seconds NOT in eg 063014 format). Suffix is removed
            self.dictionary[directory] = sortList(calList)
    
        
        #output sorted file and timeStamp lists
        '''
        self.sortList = tempList2
        self.timeStampList = tempList1
        '''
        
    def centroid(self,directory,worldCoor,refPixGuide=[629.67,318.49],refPixARCONS=[32.024,28.849]):
        '''
        Convert world coordinate into arcons coordinate.
        '''
        guidePixel = [] 
        
        fileList = timeConvert(self.dictionary[directory])
        
        for nfile in fileList:
            #include the prefix for file directory
            nfiledir = self.caldir + nfile + self.suffix
            
            #try to include distortion parameter if presented
            try:
                pix = world2pix(nfiledir,[worldCoor],'params.npz')
            except:
                pix = world2pix(nfiledir,[worldCoor], None)
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
        #eliminate any coordinate that is out of the range of the plate (0<x<44,0<y<46)
        for coor1 in arconsX:
            if coor1 > 0 and coor1 < 44 and arconsY[index] > 0 and arconsY[index] < 46:
                arconsCoor.append([coor1,arconsY[index],self.sortList[index][0:6]])
            index += 1
        
        print arconsCoor
        return arconsCoor

    def photonMapping(self,date,timeStamp,xPhotonPixel,yPhotonPixel):
        #notice here that timeStamp has to be in STRING format (eg. '065413')
        
        #directory of the files
        directory = self.caldir + date + '/'
        
        def _find_nearest(array,value):
        #find closest neighbor from a given array for a given value
            value = np.array(value)
            array = np.array(array)
            numCol = len(value)
            value = value.reshape(numCol,1)
            stackList = []
            for i in range(numCol):
                stackList.append(array)
            arraryList = np.vstack(stackList)
            #return list of index of nearest nb
            idx = (np.abs(arrayList-value)).argmin(axis=1)
            idx = idx.reshape(numCol)
            #return list of cloest nbh in seconds
            return array[idx]
            
        #Constructing lookup table if not presented
        try:
            print '> Loading lookup table in %s...' %directory
            saveName = directory + date + 'lookup.hdf5'
            h5f = h5py.File(saveName,'r')
            #loading index
            index = h5f['index'][:]
            print '> Table found in %s! Proceed to photon mapping!' %saveName
            
        except:
            raise
            print '> Lookup table not presented in %s. Constructing a new table NOW!' %directory
            
            '''
            #+0.1 is to include the very last number
            tempList = [np.arange(0,self.dimX+0.1,0.1),np.arange(0,self.dimY+0.1,0.1)]
            #pixel is the list that contains every pixels in 0.1 steps (eg, (14.3,5.3))
            pixelList = list(itertools.product(*a))
            '''
            
            #Initializing data structure. Here I use step is 0.1
            print '> Initializing structure array...'
            meshx,meshy = np.meshgrid(np.arange(0,self.dimX+0.1,0.1),np.arange(0,self.dimY+0.1,0.1))
            
            #first apply inverse transformation to each grid point. These operations are independent of the calibration file.
            deltaX = meshx - self.refPixARCONS[0]
            deltaY = -meshy + self.refPixARCONS[1]
            x = deltaX*(self.arconsPlate/self.guidePlate) + self.refPixGuide[0]
            y = deltaY*(self.arconsPlate/self.guidePlate) + self.refPixGuide[1]
            
            #here is stack in (x,y) order(but the slicing is in REVERSE order. Namely, (x=3,y=5)=arrayMatrix[5][3])
            arrayMatrix = np.dstack((x,y))
                        
            #reshape the array to pass to pix2world
            arrayYLen = arrayMatrix.shape[0]
            arrayXLen = arrayMatrix.shape[1]     
            length = arrayYLen * arrayXLen
            flatList = arrayMatrix.reshape(length,2)
            
            #initialize h5 table
            saveName = date + 'lookup.hdf5'
            h5f = h5py.File(directory+saveName,'w')
            fileList = self.dictionary[directory]
            count = 1
            for calFile in fileList:
                filePath = directory + timeConvert(calFile) + self.suffix
                print '> Calculating %s with 0.1 step (%s/%s)...' %(filePath,count,len(fileList))
                worldCoor = np.array(pix2world(filePath,flatList))
                worldCoor = worldCoor.reshape(arrayYLen,arrayXLen,2)
                #create a dataset in h5 table with its name specified by the calFile name
                h5f.create_dataset(str(calFile),data=worldCoor)
                count += 1
            
            #save the file list to h5 table for indexing purposes. Here, the file list is stored in seconds.
            h5f.create_dataset('index',data=np.array(fileList))
            
            print '> saving h5 file to %s' %directory+saveNmae
            h5f.close()
            print '> all done!'
            print '> restart the program to map the photons!'
            exit(-1)
            
                 
        
        #look up world coordinate in the table
        matchedTime = _find_nearest(index,timeConvert(timeStamp)).astype('str')
        returnList = []
        for count,listName in enumerate(matchedTime.tolist()):
            pixel = h5f[listName]
            x = xPhotonPixel[count]
            y = yPhotonPixel[count]
            RA,DEC = pixel[y][x]
            returnList.append([RA,DEC])
            
        return np.array(returnList)
        
        '''
        deltaX = xPhotonPixel - self.refPixARCONS[0]
        deltaY = -yPhotonPixel + self.refPixARCONS[1]
        
        #find (x,y) in guide pixel coordinate
        x = deltaX*(self.arconsPlate/self.guidePlate) + self.refPixGuide[0]
        y = deltaY*(self.arconsPlate/self.guidePlate) + self.refPixGuide[1]
        '''


if __name__ == '__main__':
    #testing script
    test = radec()  
    #nlist = test.centroid(worldCoor=[98.172422,-0.031578365])
    #print nlist
    #timestamp HAS TO BE IN STR FORMAT
    print test.photonMapping('20121208',['102340'],[30],[12])

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




		
