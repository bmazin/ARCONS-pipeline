import numpy as np
from sdssgaussfitter import gaussfit
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util import utils
from util.readDict import readDict
from scipy import pi
from scipy.stats import nanmedian
from mpltools	import style
#import figureHeader
import os
import sys
import glob
import tables
from util.popup import *
from photometry.LightCurve import LightCurve
from photometry.AperPhotometry import AperPhotometry

#this aperture function needs to be added to photometry section, generalized from 44x46 shape
def aperture(startpx,startpy,radius):
        r = radius
        length = 2*r 
        height = length
        allx = xrange(startpx-int(np.ceil(length/2.0)),startpx+int(np.floor(length/2.0))+1)
        ally = xrange(startpy-int(np.ceil(height/2.0)),startpy+int(np.floor(height/2.0))+1)
        mask=np.zeros((46,44))
        
        for x in allx:
            for y in ally:
                if (np.abs(x-startpx))**2+(np.abs(y-startpy))**2 <= (r)**2 and 0 <= y and y < 46 and 0 <= x and x < 44:
                    mask[y,x]=1.
        return mask

startWvl = 4000
stopWvl = 5000

scratchDir = '/home/srmeeker/scratch/PAL2014/1SWASP_J0002'

basePath = '/Scratch/DisplayStack/PAL2014/1SWASP_J0002/'
imageStackPath = os.path.join(basePath, 'ImageStacks/')

targetCentroidPath = os.path.join(basePath, 'CentroidLists/Target')
referenceCentroidPath = os.path.join(basePath, 'CentroidLists/Reference01')
targetAperturePath = os.path.join(basePath, 'ApertureStacks/Target')
referenceAperturePath = os.path.join(basePath, 'ApertureStacks/Reference01')

targetApertureRadius = 8
refApertureRadius = 7
annulusInner = 12
annulusOuter = 20

imageStacks = []
interpolatedStacks = []
times = []
targetXs = []
targetYs = []
refXs = []
refYs = []

targetCounts = []
targetSkyCounts = []

referenceCounts = []
referenceSkyCounts = []

#get centroid guess positions first
print "Grabbing Target Centroid positions..."
targetCentroidFiles = sorted(glob.glob(os.path.join(targetCentroidPath, "Centroid*.h5")))
print targetCentroidFiles
for fname in targetCentroidFiles:
    print "\n------------------------------\n"
    print fname
    centroidFile = tables.openFile(fname, mode='r')
    xPos = np.array(centroidFile.root.centroidlist.xPositions.read())
    yPos = np.array(centroidFile.root.centroidlist.yPositions.read())
    centroidFile.close()
    for i in xrange(len(xPos)):
        targetXs.append(xPos[i])
        targetYs.append(yPos[i])

#print targetXs
#print targetYs


print "Grabbing Reference Centroid positions..."
refCentroidFiles = sorted(glob.glob(os.path.join(referenceCentroidPath, "Centroid*.h5")))
for fname in refCentroidFiles:
    print "\n------------------------------\n"
    print fname
    centroidFile = tables.openFile(fname, mode='r')
    xPos = np.array(centroidFile.root.centroidlist.xPositions.read())
    yPos = np.array(centroidFile.root.centroidlist.yPositions.read())
    centroidFile.close()
    for i in xrange(len(xPos)):
        refXs.append(xPos[i])
        refYs.append(yPos[i])

#print refXs
#print refYs


print "Concatenating image stacks and JDs..."
imageFiles = sorted(glob.glob(os.path.join(imageStackPath, "ImageStack_%i-%i_flat_hp.h5"%(startWvl,stopWvl)))) #change to *_hp.npz to get hot pixel masked versions
for fname in imageFiles:
    print "\n------------------------------\n"
    print fname
    stackFile = tables.openFile(fname, mode='r')
    stackNode = stackFile.root.stack
    stack = np.array(stackNode.stack.read())
    jd = np.array(stackNode.time.read()[0])
    stackFile.close()
    print "Stack length = ", np.shape(stack)
    print "JD length = ", np.shape(jd)
    for i in xrange(len(stack[0,0,:])):
        print "Concatenating Frame ", i, " in stack"
        frame = stack[:,:,i]
        #nanMask = np.isnan(frame)
        #frame[nanMask] = 0.0
        imageStacks.append(frame)
        #try:
        #    interpFrame = utils.interpolateImage(frame, method='linear')
        #except ValueError:
        #    interpFrame = np.zeros(np.shape(frame),dtype=float)
        #    print "Empty frame encountered, filled with zeros"
        #interpolatedStacks.append(interpFrame)
        times.append(jd[i])
        #inspect some interpolated frames
        '''
        if i%20==0:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title(jd[i])
            im = ax.matshow(frame, cmap = cm.get_cmap('rainbow'), vmax = 500)
            fig.colorbar(im)
            
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title(jd[i])
            im = ax.matshow(interpFrame, cmap = cm.get_cmap('rainbow'), vmax = 500)
            fig.colorbar(im)
            plt.show()
        '''

totalAngles = 100
angleArray = np.linspace(0,2*np.pi,totalAngles)

if len(targetXs) != len(refXs):
    print "Target and Reference centroids have different number of points!!!!"
    print len(targetXs)
    print len(refXs)
    print "Exiting..."
    sys.exit(0)



#####################################################
#Aperture Phot code using AperPhotometry Class
#####################################################

targetCentroids = zip(targetXs,targetYs)
refCentroids = zip(refXs,refYs)
allCentroids = np.asarray(zip(targetCentroids,refCentroids))
expTimes = 10.0*np.ones((len(imageStacks)))
objFluxes = []
objSkyFluxes = []
refFluxes = []
refSkyFluxes = []
for i in range(len(imageStacks)):
    ap = AperPhotometry(image=imageStacks[i],centroid=allCentroids[i],expTime=expTimes,verbose=True,showPlot=True)
    fluxDict = ap.AperPhotometry(aper_radius=[9,7], sky_sub="median", annulus_inner = 12, annulus_outer = 25, interpolation="linear")
    objFluxes.append(fluxDict['flux'][0])
    refFluxes.append(fluxDict['flux'][1])
    objSkyFluxes.append(fluxDict['sky'][0])
    refSkyFluxes.append(fluxDict['sky'][1])

objFluxes = np.array(objFluxes,dtype=float)
refFluxes = np.array(refFluxes,dtype=float)
objSkyFluxes = np.array(objSkyFluxes,dtype=float)
refSkyFluxes = np.array(refSkyFluxes,dtype=float)

plt.plot(times, objFluxes, '.', color='blue')
plt.plot(times, refFluxes, '.', color='red')
plt.title("No sky sub")
plt.show()

plt.plot(times, objSkyFluxes, '.', color='blue')
plt.plot(times, refSkyFluxes, '.', color='red')
plt.title("Sky LCs")
plt.show()

plt.plot(times, objFluxes-objSkyFluxes, '.',color='blue')
plt.plot(times, refFluxes-refSkyFluxes, '.',color='red')
plt.show()

plt.plot(times, (objFluxes-objSkyFluxes)/(refFluxes-refSkyFluxes), '.',color='blue')
plt.show()


#####################################################
#PSF fitting code using Alex's double Gaussian from PSFphotometry
#####################################################
'''
verbose=True
showPlot=False
LC = LightCurve(basePath,verbose=verbose,showPlot=showPlot)
targetCentroids = zip(targetXs,targetYs)
refCentroids = zip(refXs,refYs)
allCentroids = np.asarray(zip(targetCentroids,refCentroids))
TargetLC = LC.makeLightCurve(imageStacks,allCentroids,10.0*np.ones((len(imageStacks)),dtype=int))
#RefLC=LC.makeLightCurve(imageStacks,zip(refXs,refYs),10.0*np.ones((len(images)),dtype=int))

pop(plotFunc=lambda fig,axes: axes.plot(TargetLC),title="Flux")

plt.plot(times,TargetLC[0], 'b')
plt.plot(times,TargetLC[1], 'r')
#plt.show()
plt.plot(times,TargetLC[0]/TargetLC[1], 'black')
plt.show()
'''


#####################################################
#PSF Fitting using old gaussfitter code
#####################################################
'''
paramsList = []
errorsList = []
fitImgList = []
chisqList = []
#plt.ion()

for i in xrange(len(targetXs)):
    print "\n------------------------------\n"
    print "Performing PSF fitting photometry on target frame ", i
    guessX = int(np.round(targetXs[i],0))
    guessY = int(np.round(targetYs[i],0))

    apertureMask = aperture(guessX, guessY, targetApertureRadius)
    frame = np.array(imageStacks[i])

    nanMask = np.isnan(frame)
    
    #err = np.ones(np.shape(frame))*1.0
    err = np.sqrt(frame)
    err[apertureMask==1] = np.inf #weight points closer to the expected psf higher
    frame[nanMask]=0#set to finite value that will be ignored
    err[nanMask] = np.inf#ignore these data points
    nearDeadCutoff=0#100/15 cps for 4000-6000 angstroms
    #err[frame<nearDeadCutoff] = np.inf
    entireMask = (err==np.inf)
    maFrame = np.ma.masked_array(frame,entireMask)

    guessAmp = 50.
    guessHeight = 5.
    guessWidth = 1.3
    guessParams = [guessHeight,guessAmp,guessX,guessY,guessWidth]
    limitedmin = 5*[True]
    limitedmax = 5*[True]
    minpars = [0,0,0,0,0.1] #default min pars, usually work fine
    #minpars = [0,0,27,27,1] #tighter constraint on PSF width to avoid fitting wrong peak if PSF is divided by dead pixels
    maxpars = [40,2000,43,43,10]
    usemoments=[True,True,True,True,True] #doesn't use our guess values, default
    #usemoments=[False,False,False,False,False]

    print "=========================="
    print "jd ", times[i]
    print "frame ", i
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
    print "*** outparams = ", outparams
    print "*** paramErrors = ", paramErrors

    fitimg[nanMask]=0  
    fitImgList.append(fitimg)
    frame[nanMask]=np.nan

    amp = outparams[1]
    width = outparams[4]
    xpos = outparams[2]
    ypos = outparams[3]

    flux = 2*np.pi*amp*width*width
    targetCounts.append(flux)    
    
    if (i%100==0):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(times[i])
        im = ax.matshow(imageStacks[i], cmap = cm.get_cmap('rainbow'), vmax = 500)
        ax.plot(targetXs[i]+targetApertureRadius*np.cos(angleArray), targetYs[i]+targetApertureRadius*np.sin(angleArray), 'b')
        fig.colorbar(im)
        #plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(times[i])
        im = ax.matshow(fitimg, cmap = cm.get_cmap('rainbow'), vmax = 500)
        fig.colorbar(im)
        plt.show()
    

for i in xrange(len(refXs)):
    print "\n------------------------------\n"
    print "Performing PSF fitting photometry on reference frame ", i
    guessX = int(np.round(refXs[i],0))
    guessY = int(np.round(refYs[i],0))

    apertureMask = aperture(guessX, guessY, refApertureRadius)
    frame = np.array(imageStacks[i])

    nanMask = np.isnan(frame)
    
    #err = np.ones(np.shape(frame))*1.0
    err = np.sqrt(frame)
    err[apertureMask==1] = np.inf #weight points closer to the expected psf higher
    frame[nanMask]=0#set to finite value that will be ignored
    err[nanMask] = np.inf#ignore these data points
    nearDeadCutoff=0#100/15 cps for 4000-6000 angstroms
    #err[frame<nearDeadCutoff] = np.inf
    entireMask = (err==np.inf)
    maFrame = np.ma.masked_array(frame,entireMask)

    guessAmp = 50.
    guessHeight = 5.
    guessWidth = 1.3
    guessParams = [guessHeight,guessAmp,guessX,guessY,guessWidth]
    limitedmin = 5*[True]
    limitedmax = 5*[True]
    minpars = [0,0,0,0,0.1] #default min pars, usually work fine
    #minpars = [0,0,27,27,1] #tighter constraint on PSF width to avoid fitting wrong peak if PSF is divided by dead pixels
    maxpars = [40,2000,43,43,10]
    usemoments=[True,True,True,True,True] #doesn't use our guess values, default
    #usemoments=[False,False,False,False,False]

    print "=========================="
    print "jd ", times[i]
    print "frame ", i
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
    print "*** outparams = ", outparams
    print "*** paramErrors = ", paramErrors

    fitimg[nanMask]=0  
    fitImgList.append(fitimg)
    frame[nanMask]=np.nan

    amp = outparams[1]
    width = outparams[4]
    xpos = outparams[2]
    ypos = outparams[3]

    flux = 2*np.pi*amp*width*width
    referenceCounts.append(flux)    
    
    if (i%100==0):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(times[i])
        im = ax.matshow(imageStacks[i], cmap = cm.get_cmap('rainbow'), vmax = 500)
        ax.plot(refXs[i]+refApertureRadius*np.cos(angleArray), refYs[i]+refApertureRadius*np.sin(angleArray), 'b')
        fig.colorbar(im)
        #plt.show()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(times[i])
        im = ax.matshow(fitimg, cmap = cm.get_cmap('rainbow'), vmax = 500)
        fig.colorbar(im)
        plt.show()
    

targetCounts = np.array(targetCounts, dtype=float)
#targetSkyCounts = np.array(targetSkyCounts, dtype=float)
referenceCounts = np.array(referenceCounts, dtype=float)
#referenceSkyCounts = np.array(referenceSkyCounts, dtype=float)

plt.plot(times,targetCounts, 'b')
plt.plot(times,referenceCounts, 'r')
plt.show()
plt.plot(times,(targetCounts)/(referenceCounts), 'black')
plt.show()
'''


#####################################################
#Aperture Phot code
#####################################################
'''
for i in xrange(len(targetXs)):
    print "\n------------------------------\n"
    print "Performing target aperture photometry on frame ", i
    startpx = int(np.round(targetXs[i],0))
    startpy = int(np.round(targetYs[i],0))

    apertureMask = aperture(startpx, startpy, targetApertureRadius)

    innerMask = aperture(startpx, startpy, annulusInner)
    outerMask = aperture(startpx, startpy, annulusOuter)
    annulusMask = outerMask-innerMask

    currentImage = np.array(interpolatedStacks[i])
    nanMask = np.isnan(currentImage)
    currentImage[nanMask] = 0.0#set to finite value that will be ignored

    aperturePixels = np.array(np.where(np.logical_and(apertureMask==1, nanMask==False)))
    aperturePix = aperturePixels.shape[1]
    apertureCountsPerPixel = np.sum(currentImage[aperturePixels[0],aperturePixels[1]])/aperturePix
    targetCounts.append(apertureCountsPerPixel*aperturePix)

    annulusPixels = np.array(np.where(np.logical_and(annulusMask==1, nanMask==False)))
    annulusPix = annulusPixels.shape[1]
    #annulusCountsPerPixel = np.sum(currentImage[annulusPixels[0],annulusPixels[1]])/annulusPix
    annulusCountsPerPixel = nanmedian(currentImage[annulusPixels[0],annulusPixels[1]])
    targetSkyCounts.append(annulusCountsPerPixel*aperturePix)
    
    if (i%900==0):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(times[i])
        im = ax.matshow(interpolatedStacks[i], cmap = cm.get_cmap('rainbow'), vmax = 200)
        ax.plot(targetXs[i]+targetApertureRadius*np.cos(angleArray), targetYs[i]+targetApertureRadius*np.sin(angleArray), 'b')
        ax.plot(targetXs[i]+annulusInner*np.cos(angleArray), targetYs[i]+annulusInner*np.sin(angleArray), 'r')
        ax.plot(targetXs[i]+annulusOuter*np.cos(angleArray), targetYs[i]+annulusOuter*np.sin(angleArray), 'r')
        fig.colorbar(im)
        plt.show()
    

for i in xrange(len(refXs)):
    print "\n------------------------------\n"
    print "Performing reference aperture photometry on frame ", i
    startpx = int(np.round(refXs[i],0))
    startpy = int(np.round(refYs[i],0))

    apertureMask = aperture(startpx, startpy, refApertureRadius)

    innerMask = aperture(startpx, startpy, annulusInner)
    outerMask = aperture(startpx, startpy, annulusOuter)
    annulusMask = outerMask-innerMask

    currentImage = np.array(interpolatedStacks[i])
    nanMask = np.isnan(currentImage)
    currentImage[nanMask] = 0.0#set to finite value that will be ignored

    aperturePixels = np.array(np.where(np.logical_and(apertureMask==1, nanMask==False)))
    aperturePix = aperturePixels.shape[1]
    apertureCountsPerPixel = np.sum(currentImage[aperturePixels[0],aperturePixels[1]])/aperturePix
    referenceCounts.append(apertureCountsPerPixel*aperturePix)

    annulusPixels = np.array(np.where(np.logical_and(annulusMask==1, nanMask==False)))
    annulusPix = annulusPixels.shape[1]
    #annulusCountsPerPixel = np.sum(currentImage[annulusPixels[0],annulusPixels[1]])/annulusPix
    annulusCountsPerPixel = nanmedian(currentImage[annulusPixels[0],annulusPixels[1]])
    referenceSkyCounts.append(annulusCountsPerPixel*aperturePix)
    
    if (i%900==0):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title(times[i])
        im = ax.matshow(interpolatedStacks[i], cmap = cm.get_cmap('rainbow'), vmax = 200)
        ax.plot(refXs[i]+refApertureRadius*np.cos(angleArray), refYs[i]+refApertureRadius*np.sin(angleArray), 'b')
        ax.plot(refXs[i]+annulusInner*np.cos(angleArray), refYs[i]+annulusInner*np.sin(angleArray), 'r')
        ax.plot(refXs[i]+annulusOuter*np.cos(angleArray), refYs[i]+annulusOuter*np.sin(angleArray), 'r')
        fig.colorbar(im)
        plt.show()



targetCounts = np.array(targetCounts, dtype=float)
targetSkyCounts = np.array(targetSkyCounts, dtype=float)
referenceCounts = np.array(referenceCounts, dtype=float)
referenceSkyCounts = np.array(referenceSkyCounts, dtype=float)

plt.plot(times,targetCounts,'.',color='blue')
plt.plot(times,referenceCounts,'.',color='red')
plt.title("Raw Lightcurves (no sky sub)")
plt.savefig(scratchDir+"/%i-%i_rawLCs.png"%(startWvl,stopWvl))
plt.show()

plt.plot(times, targetSkyCounts,'.',color='blue')
plt.plot(times, referenceSkyCounts,'.',color='red')
plt.title("Sky lightcurves")
plt.savefig(scratchDir+"/%i-%i_skyLCs.png"%(startWvl,stopWvl))
plt.show()

plt.plot(times,targetCounts-targetSkyCounts,'.',color='blue')
plt.plot(times,referenceCounts-referenceSkyCounts,'.',color='red')
plt.title("Sky subtracted lightcurves")
plt.savefig(scratchDir+"/%i-%i_skySubtractedLCs.png"%(startWvl,stopWvl))
plt.show()

plt.plot(times,(targetCounts)/(referenceCounts),'.',color='black')
plt.title("Target/Reference NO sky sub")
plt.savefig(scratchDir+"/%i-%i_dividedLC_NoSkySub.png"%(startWvl,stopWvl))
plt.show()

plt.plot(times,(targetCounts-targetSkyCounts)/(referenceCounts-referenceSkyCounts),'.',color='black')
plt.title("Target/Reference WITH sky sub")
plt.savefig(scratchDir+"/%i-%i_dividedLC_WithSkySub.png"%(startWvl,stopWvl))
plt.show()
'''

