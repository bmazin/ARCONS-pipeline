import pyfits
import numpy as np
from astropy import wcs
from scipy.optimize import fsolve

inputFile = 'runRecord.txt'


def _poly2(p,x1,x2):
    """
    Polynomial distortion model following SIP convention
        
    Keywords argument:
    p --- ordered coeffient list
    x1 --- x coordinate
    x2 --- y coordinate
    """
    a10,a01,a20,a02,a11,a21,a12,a30,a03 = p
    y = a10*x1 + a01*x2 + a20*x1**2 + a02*x2**2 + a11*x1*x2 + a21*(x1**2)*x2 + a12*x1*x2**2 + a30*x1**3 + a03*x2**3
    return y
    


class radec(object):

    def __init__(self,recordFile='runRecord.txt',tolError=400,caldir='./cal/'):
    
        self.sortList = []
        self.caldir = caldir
    
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
        for nfile in fileList:
            timeStamp = nfile[0:6]
            hr = int(timeStamp[0:2])
            min = int(timeStamp[2:4])
            sec = int(timeStamp[4:6])
            sec = hr*3600 + min*60 + sec
            tempList1.append(sec)
            
        tempList1.sort() 
            

        for ntime in tempList1:
            hr = ntime/3600
            min = (ntime-(hr*3600))/60
            sec = ntime-(hr*3600)-(min*60)
            
            hr = str(hr)
            min = str(min)
            sec = str(sec)
            
            if len(hr) == 1:
                hr = '0' + hr
            if len(min) == 1:   
                min = '0' + min
            if len(sec) == 1:
                sec = '0' + sec
            outputFileName = '%s%s%sfix_offCal_rotCal.fits' %(hr,min,sec)
            tempList2.append(outputFileName)        
        
        self.sortList = tempList2

    def centroid(self,worldCoor,refPixGuide=[548,370],refPixARCONS=[11.8688,12.9538],guidePlate=3.5833*10**(-5),arconsPlate=1.208*10**(-4),dimX=44,dimY=46):
        '''
        Convert world coordinate into arcons coordinate.
        '''
        
        #this function is used in FitsAnalysis also, combine them into one place
        def _world2pix(fitsImageName,world,paramFile=None):  
            #load file
            imageList = pyfits.open(fitsImageName)
            header = imageList[0].header
            #fix error
            if header['CTYPE1'] == 'RA--TAN':
                header['CTYPE1'] = 'RA---TAN'
            #coordinates conversion
            w = wcs.WCS(header)
            pix = w.wcs_world2pix(world,1)
            pix = pix.tolist()
            
            #Then if distortion params provided, find the inverse of that and apply it to find the final pixel coordinates      
            if paramFile != None:
                
                paramFile = np.load(paramFile)
                dt1 = paramFile['dt1']
                dt2 = paramFile['dt2']
                x1ref = paramFile['xref'][0]
                x2ref = paramFile['xref'][1]
                #Prevent python errors of too many files opened
                paramFile.close()
                                          
                def _equations(p,dt1,dt2,x,y,x1ref,x2ref):
                    """ 
                    Return a system of equations roots will give the reverse transformation of distortion. Remember everything is in relative coordinates.
                    """
                    x1,x2 = p
                    eqn1 = x - _poly2(dt1,x1,x2) - x1 
                    eqn2 = y - _poly2(dt2,x1,x2) - x2 
                    return (eqn1,eqn2)
                
                
                for index in range(len(pix)):
                    x = pix[index][0] - x1ref
                    y = pix[index][1] - x2ref
                
                    x1,x2 = fsolve(_equations,(x,y),args=(dt1,dt2,x,y,x1ref,x2ref))
                    
                    pix[index][0] = x1 + x1ref
                    pix[index][1] = x2 + x2ref
                    
            return pix

        guidePixel = [] 

        for nfile in self.sortList:
            nfiledir = self.caldir + nfile
            try:
                pix = _world2pix(nfiledir,[worldCoor],'params.npz')
            except:
                pix = _world2pix(nfiledir,[worldCoor], None)
            guidePixel.append([nfile[0:6],pix[0]]) 

       
        x = np.array([coorX for coorX in (guidePixel[index][1][0] for index in range(len(guidePixel)))])

        y = np.array([coorY for coorY in (guidePixel[index][1][1] for index in range(len(guidePixel)))])

        deltaX = (x-refPixGuide[0])*guidePlate/arconsPlate
        deltaY = (y-refPixGuide[1])*guidePlate/arconsPlate

        arconsX = refPixARCONS[0] + deltaX
        arconsY = refPixARCONS[1] - deltaY

        index = 0
        arconsCoor = []
        print arconsX
        for coor1 in arconsX:
            if coor1 > 0 and coor1 < 44 and arconsY[index] > 0 and arconsY[index] < 46:
                arconsCoor.append([coor1,arconsY[index],self.sortList[index][0:6]])
            index += 1
        return arconsCoor

    def photonMapping(self):
        pass


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




		
