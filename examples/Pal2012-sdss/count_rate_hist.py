from sdssgaussfitter import gaussfit
import numpy as np
from util import utils
import matplotlib.pyplot as plt

def aperture(startpx,startpy,radius=3):
    r = radius
    length = 2*r 
    height = length
    allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
    ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
    pixx = []
    pixy = []
    mask=np.ones((46,44))
    for x in allx:
        for y in ally:
            if (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 <= (r)**2 and 0 <= y and y < 46 and 0 <= x and x < 44:
                mask[y,x]=0.
    return mask

def gaussian(height, center_x, center_y, width_x, width_y,offset):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)+offset

#testy = np.array([[gaussian(2,10,10,3,3,5)(x,y) for y in range(46)] for x in range(44)])
#utils.plotArray(testy,cbar=True)
stackDict = np.load('/home/pszypryt/sdss_data/20121208/RedTest-ImageStack.npz')
stack = stackDict['stack']
paramsList = []
fitImgList = []
frames = []
for iFrame in range(0,np.shape(stack)[2]):
    frame = stack[:,:,iFrame]
    plt.hist(np.ravel(frame),bins=300,range=(0,1000))
    plt.show()

    nanMask = np.isnan(frame)
    frame[nanMask] = 0
    #frame = np.ma.masked_array(frame,mask=nanMask)
    #frames.append(frame)
    #utils.plotArray(frame,cbar=True)

    
#utils.makeMovie(frames,cbar=True,outName='nlttImageStackNorm.gif',normMax=5000)
