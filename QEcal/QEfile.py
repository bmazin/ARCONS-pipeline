'''
Author: Alex Walter
Date: 5-13-2013

Helper class for QECalibration.py

We do the actual QE calculation here
'''


from os.path import isfile
import numpy as np
import math

class QEfile():

    def __init__(self, fn=None, ang=0):
        if not isfile(fn):
            raise fileNameError(fn)
        self.filename=fn
        self.angle=ang
        self.loadData()
        self.areaArcons = (222.0*10**-6)**2             #[m^2] effective area of arcons
                                                        #area of inductor * magnification of microlense = area of pixel (approximately)
        self.areaDect = 1.0*10.0**-4.                    #[m^2] area of optical detector
        self.areaIRDect = 1.0*math.pi*(1.5*10**-3)**2   #[m^2] area of IR detector (wavelength >= 1100
        self.magnification = 1.2                        # magnification with Palamar optics
        

    def findQE(self, nArcons, wavelength):
        #nArcons: average number of photons detected at arcons per second
        #nDect: average number of photons detected at power meter per second
        #wavelength in nanometers
       
        nDect = self.data[np.where(self.data[:,0]==wavelength)[0][0],-7]
        nDect*=10.0**7    #[photons/s]
#        print 'Detected Count Rate: ' + str(nDect)
#        print 'Arcons Count Rate: ' + str(nArcons)

        if wavelength >= 1099.999:
            dectArea=self.areaIRDect
        else:
            dectArea=self.areaDect
        
#        print 'Detector Area: ' + str(dectArea)
        return 1.0*nArcons/(self.magnification**2*self.areaArcons)/(nDect/(dectArea))

    def loadData(self):
        print 'loading data from '+self.filename
        d=np.loadtxt(str(self.filename))
        self.data=np.zeros((len(d),len(d[0])+6))
        self.data[:,:-6]=d
        self.determineTiming()

    def saveParam(self):
        fileparam = 'param'.join(str(self.filename).rsplit('txt',1))
        np.savetxt(fileparam, self.data[:,-6:],fmt='%i', delimiter='\t')
        print 'saved parameters: ' + str(fileparam)

    def determineTiming(self):
        #Can make a more sophisticated guess later
        #Works well for 20130503POL2.txt, obs_20130503-234952.h5

        #load data from param file
        fileparam = 'param'.join(str(self.filename).rsplit('txt',1))
        if isfile(fileparam):
            param = np.loadtxt(str(fileparam))
            self.data[:,-6:] = param
            print 'loaded parameters: ' + str(fileparam)
        else:
            startTime=125
            widthTime=15
            troughTime=70
            start2Num=15
            start2Time=1465

            for i in range(len(self.data)):
                f1=startTime
                j=i
                if i>=(start2Num-1):
                    f1=start2Time
                    j=i-(start2Num-1)
                
                f1+=j*(widthTime+troughTime)
                f2=f1+widthTime
                t1=(3*f1-f2-troughTime)/2.0
                t2=(f2+f1-troughTime)/2.0
                t3=(f1+f2+troughTime)/2.0
                t4=(3*f2-f1+troughTime)/2.0
                self.data[i,-6]=round(t1)
                self.data[i,-5]=round(t2)           
                self.data[i,-4]=round(f1)
                self.data[i,-3]=round(f2)
                self.data[i,-2]=round(t3)
                self.data[i,-1]=round(t4)
            self.saveParam()

#    def getData(self):
#        return self.data

#    def setData(self, x=-1, y=-1, value=0):
#        print 'old value: ' +str(self.data[x,y])
#        print 'new value: '+str(value)
#        self.data[x,y] = value

#    def getLength(self):
#        return len(self.data)

#    def getAngle(self):
#        return self.angle
        
        



class fileNameError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value) + ' is not a file'




