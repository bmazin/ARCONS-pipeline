from sdssgaussfitter import gaussfit
import numpy as np
from util import utils
import matplotlib.pyplot as plt
import sys

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

stackDict = np.load('/Scratch/dataProcessing/SDSS_J0926/AllData/Dec11SIImageStackRed.npz')
stack = stackDict['stack']
if len(sys.argv) == 1:
    print 'Useage: ',sys.argv[0],' iFrame'
    print """
    set0 Frames 0-179
    """
    exit(1)
iFrame = int(sys.argv[1])
frame = stack[:,:,iFrame]
#    plt.hist(np.ravel(frame),bins=100,range=(0,5000))
#    plt.show()

nanMask = np.isnan(frame)
frame[nanMask] = 0
frame = np.ma.masked_array(frame,mask=nanMask)
utils.plotArray(frame,cbar=True)

    
