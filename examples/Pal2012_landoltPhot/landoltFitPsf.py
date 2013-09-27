from sdssgaussfitter import gaussfit
import numpy as np
from scipy import interpolate
import os,sys
from util import utils
from util.readDict import readDict
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def aperture(startpx,startpy,radius=7):
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

#def gaussian(height, center_x, center_y, width_x, width_y,offset):
#    """Returns a gaussian function with the given parameters"""
#    width_x = float(width_x)
#    width_y = float(width_y)
#    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)+offset

#testy = np.array([[gaussian(2,10,10,3,3,5)(x,y) for y in range(46)] for x in range(44)])
#utils.plotArray(testy,cbar=True)

param = readDict()
#param.read_from_file('G158-100params.dict')
#param.read_from_file('pg0220params.dict')
#param.read_from_file('landolt9542params.dict')
#param.read_from_file('corot18params.dict')
if len(sys.argv)<2:
    print "Provide file name to fit.  Syntax >>python fitPsf.py objectparams.dict [filenumber]"
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
#NumFiles = param['NumFiles']
#for filenum in range(len(NumFiles)):
#    if NumFiles[filenum] > 0:
#        NumFiles[filenum] = NumFiles[filenum]*FramesPerFile
#NumFrames = NumFiles

NumFrames = 31
print NumFrames
guessX = param['guessX'][0]
guessY = param['guessY'][0]

stackDict = np.load(npzLoadFile)

stack = stackDict['stack']
wvls = stackDict['wvls']
print len(wvls)
paramsList = []
errorsList = []
fitImgList = []
chisqList = []
plt.ion()

for iFrame in range(0,np.shape(stack)[0]):
    frame = stack[iFrame,:,:]

    #for interval in xrange(len(NumFrames)-1):
    #    if NumFrames[interval] != NumFrames[interval+1]:
    #        if NumFrames[interval] < iFrame <= NumFrames[interval+1]:
    #            guessX = guessX[interval]
    #            guessY = guessY[interval]
    #             print guessX, guessY
 
    #TRY SPLINE INTERPOLATION TO FILL IN BLANK SPACES AND HELP FIT
    grid = np.zeros((np.shape(frame)[0],np.shape(frame)[1],2),dtype = int)
    for i in xrange(np.shape(frame)[0]):
        for j in xrange(np.shape(frame)[1]):
            grid[i,j]= i,j
    #print interpFrame
    reshapeFrame =np.reshape(frame,(46*44))
    reshapeGrid = np.reshape(grid,(46*44,2))
    interpFrame = interpolate.griddata(reshapeGrid[reshapeFrame!=0],reshapeFrame[reshapeFrame!=0],reshapeGrid)
    print np.reshape(grid,(46*44,2))
    interpFrame = np.reshape(interpFrame,(np.shape(frame)[0],np.shape(frame)[1]))
    print frame
    origFrame = frame
    print interpFrame

    frame = interpFrame

    nanMask = np.isnan(frame)

    apertureMask = aperture(guessX,guessY,radius=7)
    #err = np.sqrt(frame) #divide by 2 to constrain PSF fit even tighter to avoid fitting to wrong peak if PSF is divided by dead pixels
    err = np.ones(np.shape(frame))
    err[frame<10] = 100
    frame[nanMask]=0#set to finite value that will be ignored
    #err[nanMask] = 1E6 #ignore these data points
    err[frame==0] = 1E6
    #err[apertureMask==1] = 1 #np.sqrt(frame[apertureMask==1]) #weight points closer to the expected psf higher
    
    nearDeadCutoff=1#100/15 cps for 4000-6000 angstroms
    err[frame<nearDeadCutoff] = 1E6
    entireMask = (err==1E6)
    maFrame = np.ma.masked_array(frame,entireMask)
    guessAmp = 120.
    guessHeight = 3.
    guessWidth=1.6
    guessParams = [guessHeight,guessAmp,guessX,guessY,guessWidth]
    limitedmin = 5*[True] 
    limitedmax = 5*[True]
    #minpars = [0,0,0,0,.1] #default min pars, usually work fine
    minpars = [0,116,25,25,1] #tighter constraint on PSF width to avoid fitting wrong peak if PSF is divided by dead pixels
    maxpars = [3,200,40,40,10]
    #usemoments=[True,True,True,True,True] #doesn't use our guess values, default
    usemoments=[False,False,False,False,False]
    
    out = gaussfit(data=maFrame,err=err,params=guessParams,returnfitimage=True,quiet=True,limitedmin=limitedmin,limitedmax=limitedmax,minpars=minpars,maxpars=maxpars,circle=1,usemoments=usemoments,returnmp=True)
    mp = out[0]

    outparams = mp.params
    paramErrors = mp.perror
    chisq = mp.fnorm
    dof = mp.dof
    reducedChisq = chisq/dof
    print "reducedChisq =", reducedChisq
    fitimg = out[1]
    chisqList.append([chisq,dof])


    paramsList.append(outparams)
    errorsList.append(paramErrors)
    print "outparams = ", outparams
    print "paramErrors = ", paramErrors

#    expectedResiduals = np.ma.masked_array(np.sqrt(frame),mask=entireMask)
#    residuals = np.ma.masked_array(np.abs(frame-fitimg),mask=entireMask)
#    utils.plotArray(expectedResiduals,cbar=True)
#    utils.plotArray(residuals,cbar=True)
#    fig = plt.figure()
#    ax = fig.add_subplot(111,projection='3d')
#    x = np.arange(0,44)
#    y = np.arange(0,46)
#    X,Y = np.meshgrid(x,y)
#    linearMask = np.ravel(entireMask==0)
#    ax.plot_wireframe(X,Y,fitimg)
#    ax.scatter(outparams[2],outparams[3],outparams[0]+outparams[1],c='black')
#    ax.scatter(np.ravel(X)[linearMask],np.ravel(Y)[linearMask],np.ravel(frame)[linearMask],c='red')
#
    fitimg[nanMask]=0
#    print fitimg[np.isnan(fitimg)]            
    fitImgList.append(fitimg)

#    utils.plotArray(frame,cbar=True)
#    utils.plotArray(maFrame,cbar=True)
#    utils.plotArray(fitimg,cbar=True)
#    plt.show()
#    utils.confirm('Enter to continue.')
#    plt.close()
#    plt.close()
#    plt.close()

    frame[nanMask]=np.nan


#    fig = plt.figure()
#    ax1=fig.add_subplot(211)
#    ax2 = fig.add_subplot(212)
#    for iRow in range(len(frame)):
#        ax1.scatter(range(44),frame[iRow,:],c='red',marker='o',alpha=.5,label='data')
#        ax1.scatter(range(44),fitimg[iRow,:],c='blue',marker='^',alpha=.5,label='fit')
#    ax1.set_title('Fit seen along Cols')
#    for iCol in range(np.shape(frame)[1]):
#        ax2.scatter(range(46),frame[:,iCol],c='red',marker='o',alpha=.5,label='data')
#        ax2.scatter(range(46),fitimg[:,iCol],c='blue',marker='^',alpha=.5,label='fit')
#    ax2.set_title('Fit seen along Rows')
#    plt.show()

plt.close()
print 'closed'
cube = np.array(fitImgList)
chisqs = np.array(chisqList)
params = np.array(paramsList)
errors = np.array(errorsList)

np.savez(npzfitpsf,fitImg=cube,params=params,errors=errors,chisqs=chisqs,wvls=wvls)
print 'saved'
#utils.makeMovie(fitImgList,frameTitles=wvls, cbar=True, outName=giffitpsf, normMin=0, normMax=50)
stop = 'n'
while stop != 'q':
    plt.matshow(interpFrame,vmin=0,vmax = 100)
    plt.matshow(origFrame,vmin=0,vmax=100)
    plt.show()
    stop = raw_input(" q for stop ----> ")



