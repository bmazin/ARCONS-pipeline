from sdssgaussfitter import gaussfit
import numpy as np
from util import utils
from mpl_toolkits.mplot3d import Axes3D
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

#def gaussian(height, center_x, center_y, width_x, width_y,offset):
#    """Returns a gaussian function with the given parameters"""
#    width_x = float(width_x)
#    width_y = float(width_y)
#    return lambda x,y: height*np.exp(-(((center_x-x)/width_x)**2+((center_y-y)/width_y)**2)/2)+offset

#testy = np.array([[gaussian(2,10,10,3,3,5)(x,y) for y in range(46)] for x in range(44)])
#utils.plotArray(testy,cbar=True)

stackDict = np.load('/Scratch/dataProcessing/SDSS_J0926/AllData/AllDataSIImageStackRednewCal.npz')

stack = stackDict['stack']
jd = stackDict['jd']
paramsList = []
errorsList = []
fitImgList = []
chisqList = []
plt.ion()
for iFrame in range(0,np.shape(stack)[2]):
    frame = stack[:,:,iFrame]
    nanMask = np.isnan(frame)
    if 1700 < iFrame < 3100:
        guessX = 31
        guessY = 30
    if iFrame > 4800:
        guessX = 16
        guessY = 15
    else:
        guessX = 14
        guessY = 8

    apertureMask = aperture(guessX,guessY,radius=10)
    err = np.sqrt(frame)
    #err = np.ones(np.shape(frame))
    err[apertureMask==1] = np.inf#weight points closer to the expected psf higher
    frame[nanMask]=0#set to finite value that will be ignored
    err[nanMask] = np.inf#ignore these data points
    nearDeadCutoff=1#100/15 cps for 4000-6000 angstroms
    err[frame<nearDeadCutoff] = np.inf
    entireMask = (err==np.inf)
    maFrame = np.ma.masked_array(frame,entireMask)
    guessAmp = 600.
    guessHeight = 675.
    guessWidth=3.
    guessParams = [guessHeight,guessAmp,guessX,guessY,guessWidth]
    limitedmin = 5*[True] 
    limitedmax = 5*[True]
    minpars = [0,0,0,0,.1]
    maxpars = [5000,10000,43,45,10]
    usemoments=[True,True,True,True,True] #doesn't use our guess values
    
    
    out = gaussfit(data=maFrame,err=err,params=guessParams,returnfitimage=True,quiet=True,limitedmin=limitedmin,limitedmax=limitedmax,minpars=minpars,maxpars=maxpars,circle=1,usemoments=usemoments,returnmp=True)
    mp = out[0]

    outparams = mp.params
    paramErrors = mp.perror
    chisq = mp.fnorm
    dof = mp.dof
    reducedChisq = chisq/dof
    print reducedChisq
    fitimg = out[1]
    chisqList.append([chisq,dof])


    paramsList.append(outparams)
    errorsList.append(paramErrors)
    print outparams,paramErrors

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
    print fitimg[np.isnan(fitimg)]            
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

np.savez('/Scratch/dataProcessing/SDSS_J0926/AllData/AllDataSIfitpsfRednewCal.npz',fitImg=cube,params=params,errors=errors,chisqs=chisqs,jd=jd)
print 'saved'
utils.makeMovie(fitImgList,cbar=True,outName='/Scratch/dataProcessing/SDSS_J0926/AllData/AllDataSIfitpsfRednewCal.gif')
#,normMax=1000)

