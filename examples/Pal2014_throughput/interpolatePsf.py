from sdssgaussfitter import gaussfit
import numpy as np
import os,sys
from util import utils
from util.readDict import readDict
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import interpolate
from scipy.optimize.minpack import curve_fit
from numpy import exp


param = readDict()
if len(sys.argv)<2:
    print "Provide file name to fit.  Syntax >>python interpolatePsf.py objectparams.dict [filenumber]"
    sys.exit(1)

#read in parameter file as command line argument
param.read_from_file(sys.argv[1])

#provide optional file number if the object in the param file has alternate .npz files to be specified individually
fileNum = None
if len(sys.argv)>2:
    fileNum = "_"+str(sys.argv[2])

npzLoadFile = param['npzLoadFile']
npzfitpsf = param['npzfitpsf']
giffitpsf = param['giffitpsf']

if fileNum != None:
    npzLoadFile = npzLoadFile.split('.')[0]+fileNum+'.'+npzLoadFile.split('.')[1]
    npzfitpsf = npzfitpsf.split('.')[0]+fileNum+'.'+npzfitpsf.split('.')[1]
    giffitpsf = giffitpsf.split('.')[0]+fileNum+'.'+giffitpsf.split('.')[1]

FramesPerFile = param['FramesPerFile']

stackDict = np.load(npzLoadFile)
stack = stackDict['stack']
wvls = stackDict['wvls']
fitImgList = []

for iFrame in range(0,np.shape(stack)[0]):
    frame = stack[iFrame,:,:]
    finalshape = np.shape(frame)

    frame[frame>70] = 70

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(wvls[iFrame])
    im = ax.matshow(frame, cmap = cm.get_cmap('rainbow'))
    fig.colorbar(im)
    #plt.show()

    dataPoints = np.where(frame!=0) #data points are only pixels with counts
    D,n = np.shape(dataPoints)
    dataPoints = np.array((dataPoints[0],dataPoints[1]),dtype=np.int).transpose() #put into format griddata wants

    interpPoints = np.where(frame!=np.nan) #final image is full 44x46 array
    D,n = np.shape(interpPoints)
    interpPoints = np.array((interpPoints[0],interpPoints[1]),dtype=np.int).transpose() #put into format griddata wants

    intFrame = interpolate.griddata(dataPoints, frame[np.where(frame!=0)], interpPoints, 'linear')
    intFrame = np.reshape(intFrame, finalshape)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(wvls[iFrame])
    im = ax.matshow(intFrame, cmap = cm.get_cmap('rainbow'))#, vmax = 40)
    fig.colorbar(im)
    plt.show()

    fitImgList.append(intFrame)

cube = np.array(fitImgList)
np.savez(npzfitpsf,cube=cube,wvls=wvls)
print 'saved'
utils.makeMovie(fitImgList,frameTitles=wvls, cbar=True, outName=giffitpsf, normMin=0, normMax=50)


