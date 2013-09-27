from sdssgaussfitter import gaussfit
import numpy as np
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
print "There should be this many frames: ", NumFrames
guessX = param['guessX'][0]
guessY = param['guessY'][0]

stackDict = np.load(npzLoadFile)

stack = stackDict['stack']
wvls = stackDict['wvls']
print "The file actually has this many: ", len(wvls)
paramsList = []
errorsList = []
fitImgList = []
chisqList = []
plt.ion()



for iFrame in range(0,np.shape(stack)[0]):
    frame = stack[iFrame,:,:]
    
    #print "Frame max= ", np.nanmax(frame,axis=None)
    #frame *= CorrFactors
    #print "Corrected Frame max= ", np.nanmax(frame,axis=None)

    nanMask = np.isnan(frame)

    #for interval in xrange(len(NumFrames)-1):
    #    if NumFrames[interval] != NumFrames[interval+1]:
    #        if NumFrames[interval] < iFrame <= NumFrames[interval+1]:
    #            guessX = guessX[interval]
    #            guessY = guessY[interval]
    #             print guessX, guessY
    '''
    apertureMask = aperture(guessX,guessY,radius=4)
    err = np.sqrt(frame) #divide by 2 to constrain PSF fit even tighter to avoid fitting to wrong peak if PSF is divided by dead pixels
    err[frame > 100]=np.inf
    #err[frame<10] = 100
    frame[nanMask] = 0 #set to finite value that will be ignored
    err[nanMask] = np.inf #ignore these data points
    err[frame==0] = np.inf
    err[apertureMask==1] = 1.0 #np.sqrt(frame[apertureMask==1])/2.0 #weight points closer to the expected psf higher
    
    nearDeadCutoff = 1 #100/15 cps for 4000-6000 angstroms
    err[frame<nearDeadCutoff] = np.inf
    entireMask = (err==np.inf)
    maFrame = np.ma.masked_array(frame,entireMask)
    '''
    apertureMask = aperture(guessX,guessY,radius=7)
    
    #if iFrame < 19:
    #    err = np.ones(np.shape(frame))*10.0
    #else:
    #    err = np.zeros(np.shape(frame))
    err = np.ones(np.shape(frame))*10.0
    err[apertureMask==1] = np.inf #weight points closer to the expected psf higher
    err[frame>100] = np.inf
    frame[nanMask]=0#set to finite value that will be ignored
    err[nanMask] = np.inf#ignore these data points
    nearDeadCutoff=1#100/15 cps for 4000-6000 angstroms
    err[frame<nearDeadCutoff] = np.inf
    entireMask = (err==np.inf)
    maFrame = np.ma.masked_array(frame,entireMask)

    guessAmp = 30.
    guessHeight = 5.
    guessWidth = 1.3
    guessParams = [guessHeight,guessAmp,guessX,guessY,guessWidth]
    limitedmin = 5*[True]
    limitedmax = 5*[True]
    minpars = [0,0,0,0,0.1] #default min pars, usually work fine
    #minpars = [0,0,27,27,1] #tighter constraint on PSF width to avoid fitting wrong peak if PSF is divided by dead pixels
    maxpars = [40,200,43,43,10]
    #maxpars = [40,200,33,33,10]
    if iFrame == 27:
        minpars = [8,5,0,0,0.5]
        maxpars = [30,25,43,45,1.1]
    if iFrame == 28:
        minpars = [8,5,0,0,0.5]
        maxpars = [30,25,43,45,1.1]
    if iFrame == 29:
        minpars = [8,5,0,0,0.5]
        maxpars = [30,25,43,45,1.1]
    if iFrame == 30:
        minpars = [8,5,0,0,0.5]
        maxpars = [30,25,43,45,1.10]
    
    usemoments=[True,True,True,True,True] #doesn't use our guess values, default
    #usemoments=[False,False,False,False,False]

    print "=========================="
    print wvls[iFrame]
    print "frame ",iFrame
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
utils.makeMovie(fitImgList,frameTitles=wvls, cbar=True, outName=giffitpsf, normMin=0, normMax=50)




















