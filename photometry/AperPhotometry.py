'''
Author: Seth Meeker
Date: Jan 7, 2014

This code performs Aperture photometry on a target and any number of guide stars in the field
'''

import warnings
import tables
from tables import *
import os
from scipy.stats import nanmedian
import itertools

from photometry.Photometry import Photometry
from photometry.plot3DImage import *
from util.ObsFile import ObsFile
from util import utils
from util.FileName import FileName
from util.fitFunctions import model_list
from util.mpfit import mpfit
from util.popup import *

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
#from http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

#from http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


class AperPhotometry(Photometry):

    def __init__(self,image,centroid,expTime=None,verbose=False,showPlot=False):
        '''
        Inputs:
            image - 2D array of data (0 for dead pixel, shouldn't be any nan's or infs)
                  - Should be fully calibrated, dead time corrected, and scaled up to the effective integration time
            centroid - list of (col,row) tuples. The first tuple is the target location. The next are reference stars in the field 
            expTime - 2d array of effective exposure times (same size as image)
            verbose - show error messages
            showPlot - show and pause after each frame
        '''
        self.verbose=verbose
        self.showPlot=showPlot
        
        super(AperPhotometry,self).__init__(image=image,centroid=centroid,expTime=expTime)

    def AperPhotometry(self, aper_radius=5, sky_sub="median", annulus_inner=10, annulus_outer=15, interpolation="linear"):
        '''
        Inputs:
            aper_radius - double or list of doubles of the same length as self.centroid. Number of pixels around the star to be used in aperture
            sky_sub - indicates type of sky subtraction to be performed.
                      "median" - to use median sky value in annulus
                      "fit" - to mask objects and estimate sky in aperture with polynomial fit
            annulus_inner - double or list of doubles of same length as self.centroid. Gives radius of inner part of sky annulus for use in "median" sky sub. Ignored if sky_sub is not "median" [target_apertureRad, ref0_apertureRad, ...]
            annulus_outer - double or list of doubles of same length as self.centroid. Gives radius of outer part of sky annulus for use in "median" sky sub. Ignored if sky_sub is not "median"
            interpolation - select type of interpolation to perform on image before doing photometry
                            None - to skip interpolation
                            "cubic" - to do 2D cubic interpolation
                            "linear" - to do 2D linear interpolation (default)
                            
        Returns: Dictionary with keywords. These keywords should be the same as in aperPhotometryDataDescription in headers.DisplayStackHeaders
            flux - array of flux values. [target_flux, target_sky_flux, ref0_flux, ...]
            skyFlux - array of sky flux values, scaled to same n_pix as object flux. [target_sky, ref0_sky, ...]
            flag - flag to indicate successful analysis. 0 is good
            apertureRad - same as input
            annulusInnerRad - same as input
            annulusOuterRad - same as input
            interpolation - same as input
            
            
        '''

        flux = []
        sky = []
        apRadOut = []
        annInOut = []
        annOutOut = []
        flag = 0

        #if sky fitting is selected for sky subtraction, do masking and produce polynomial sky image
        if sky_sub == "fit":
            warnings.warn("Sky fitting not well debugged. Use with extreme caution",UserWarning)
            skyIm = self.fitMaskedSky(self.image, aper_radius)

        if interpolation != None:
            if self.verbose:
                print "Performing %s interpolation on image"%interpolation
            try:
                interpImage = utils.interpolateImage(self.image, method=interpolation)
            except ValueError:
                interpImage = np.zeros(np.shape(self.image),dtype=float)
                print "Empty frame encountered on interpolation, filled with zeros"
            if self.showPlot:
                plotArray(title='Before interpolation', image=self.image)
                plotArray(title='After interpolation', image=interpImage)
        else:
            interpImage = self.image

        if self.showPlot:
            totalAngles = 100
            angleArray = np.linspace(0,2*np.pi,totalAngles)
            form=PopUp(title='Aperture photometry',showMe=False)
            form.plotArray(interpImage)

        #step through each star in centroid list: [target, ref0, ref1, ...]
        for star_i in range(len(self.centroid)):
                try: #check if different aperture radii are set for each star
                    radius=aper_radius[star_i]
                except TypeError:
                    radius=aper_radius

                try: #check if different annulus radii are set for each star
                    ann_in=annulus_inner[star_i]
                except TypeError:
                    ann_in=annulus_inner

                try: #check if different annulus radii are set for each star
                    ann_out=annulus_outer[star_i]
                except TypeError:
                    ann_out=annulus_outer

                objectFlux, nObjPix = self.getApertureCounts(interpImage,radius,self.centroid[star_i])
                if self.showPlot: #add aperture to output image if showPlot is true
                    form.axes.plot(self.centroid[star_i][0]+radius*np.cos(angleArray), self.centroid[star_i][1] + radius*np.sin(angleArray), 'b')

                if sky_sub == "fit":
                    skyFlux, nSkyPix = self.getApertureCounts(skyIm,radius,self.centroid[star_i])
                    ann_in = 0
                    ann_out = 0
                else:
                    skyFlux, nSkyPix = self.getAnnulusCounts(interpImage,ann_in,ann_out,self.centroid[star_i])
                    if self.showPlot: #add annulus to output image if showPlot is true
                        form.axes.plot(self.centroid[star_i][0]+ann_in*np.cos(angleArray), self.centroid[star_i][1] + ann_in*np.sin(angleArray), 'r')
                        form.axes.plot(self.centroid[star_i][0]+ann_out*np.cos(angleArray), self.centroid[star_i][1] + ann_out*np.sin(angleArray), 'r')
                flux.append(objectFlux)
                sky.append(skyFlux*(float(nObjPix)/float(nSkyPix)))

                #output radii used for aperture and annuli for bookkeeping
                apRadOut.append(radius)
                annInOut.append(ann_in)
                annOutOut.append(ann_out)

        if self.showPlot:
            form.show()

        return {'flux': np.asarray(flux), 'skyFlux': np.asarray(sky), 'apertureRad':np.asarray(apRadOut), 'annulusInnerRad':np.asarray(annInOut), 'annulusOuterRad':np.asarray(annOutOut), 'flag':flag, 'interpolation':interpolation}


    def getApertureCounts(self, im, radius, center):
        startpx = int(np.round(center[0]))
        startpy = int(np.round(center[1]))
        apertureMask = aperture(startpx, startpy, radius)    
        nanMask = np.isnan(im)
        im[nanMask] = 0.0#set to finite value that will be ignored
        aperturePixels = np.array(np.where(np.logical_and(apertureMask==1, nanMask==False)))
        nApPix = aperturePixels.shape[1]
        apertureCounts = np.sum(im[aperturePixels[0],aperturePixels[1]])
        if self.verbose:
            print "Aperture Counts = ", apertureCounts
            print "Aperture pixels = ", nApPix
        return [apertureCounts,nApPix]

    def getAnnulusCounts(self, im, annulusInner, annulusOuter, center):
        startpx = int(np.round(center[0]))
        startpy = int(np.round(center[1]))
        innerMask = aperture(startpx, startpy, annulusInner)
        outerMask = aperture(startpx, startpy, annulusOuter)
        annulusMask = outerMask-innerMask
        nanMask = np.isnan(im)
        annulusPixels =  np.array(np.where(np.logical_and(annulusMask==1, nanMask==False)))
        nAnnPix = annulusPixels.shape[1]
        annulusCounts = nanmedian(im[annulusPixels[0],annulusPixels[1]])*nAnnPix
        if self.verbose:
            print "Annulus Counts = ", annulusCounts
            print "Annulus pixels = ", nAnnPix
        return [annulusCounts, nAnnPix]

    def fitMaskedSky(self, image, aper_radius):
        '''
        fit polynomial to sky background to estimate sky under aperture
        '''
        im = np.array(image)
        nRows,nCols = np.shape(im)
        nanMask = np.isnan(im)
        xx,yy = np.meshgrid(np.arange(nCols),np.arange(nRows))

        #mask out object apertures
        apertureMask=np.zeros((np.shape(im)),dtype=int)
        for i in range(len(self.centroid)):
            center=self.centroid[i]
            try: #check if different aperture radii are set for each star
                radius=aper_radius[i]
            except TypeError:
                radius=aper_radius
            startpx = int(np.round(center[0]))
            startpy = int(np.round(center[1]))
            apertureMask += aperture(startpx, startpy, radius)
        
        maskedPix = np.array(np.where(apertureMask>0))
        im[maskedPix[0], maskedPix[1]]=0

        if self.showPlot:
            plotArray(title='Masked sky', image=im)

        z = im.ravel()
        x = xx.ravel()
        y = yy.ravel()
        x = x[z != 0]
        y = y[z != 0]
        z = z[z != 0]

        try:
            PolyFitCoeffs = polyfit2d(x,y,z)
            # Evaluate it on a grid...
            fitIm = polyval2d(xx, yy, PolyFitCoeffs)
        except ValueError:
            fitIm = np.zeros(np.shape(im),dtype=float)
            print "Empty frame encountered on interpolation, filled with zeros"

        if self.showPlot:
            plotArray(title='Fitted sky Image', image=fitIm)

        return fitIm









