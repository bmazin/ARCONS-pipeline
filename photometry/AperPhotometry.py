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

#needs to be updated
def readPhotometryFile(filename):
    photometryFile = tables.openFile(filename, mode='r')
    photometryTable = photometryFile.getNode(photometryFile.root,'AperPhotometry')._f_getChild('photometry')
    startTimes = photometryTable.col('time')
    intTimes = photometryTable.col('integration_time')
    flux = photometryTable.col('flux')
    sky = photometryTable.col('sky')
    apertureRad = photometryTable.col('apertureRad')
    annulusInnerRad = photometryTable.col('annulusInnerRad')
    annulusOuterRad = photometryTable.col('annulusOuterRad')
    flag = photometryTable.col('flag')
    photometryFile.close()
    
    return {'startTimes': startTimes, 'intTimes': intTimes, 'flux': flux, 'sky': sky, 'apertureRad': apertureRad, 'annulusInnerRad': annulusInnerRad, 'annulusOuterRad': annulusOuterRad, 'flag': flag}


def writePhotometryFile(fluxDict_list, im_dict, filename,verbose = False):
    if os.path.isfile(filename):
        warnings.warn("Overwriting photometry file: "+str(filename),UserWarning)
        
    photometryFile = tables.openFile(filename, mode='w')
    photometryGroup = photometryFile.createGroup(photometryFile.root, 'AperPhotometry', 'Table of aperture photometry data')

    num_stars = int(len(fluxDict_list[0]['flux']))
    AperPhotometry_Description = {
        "time"              : Float64Col(),             # Start time of images
        "integration_time"  : Float64Col(),             # integration time of individual images
        "flux"              : Float64Col(num_stars),    # array of object flux values. [target_flux, ref0_flux, ...]
        "sky"               : Float64Col(num_stars),    # array of sky flux values, scaled to same n_pix as object flux. [target_sky, ref0_flux, ...]
        "apertureRad"       : Float64Col(num_stars),    # radius of aperture used for each object [target_apertureRad, ref0_apertureRad, ...]
        "annulusInnerRad"   : Float64Col(num_stars),    # radius of inner annulus used for sky subtraction. 0s if sky fitting was used.
        "annulusOuterRad"   : Float64Col(num_stars),    # radius of outer annulus used for sky subtraction. 0s if sky fitting was used
        "flag"              : UInt16Col()}              # flag to indicate if sky fit is good (0) 
        
    photometryTable = photometryFile.createTable(photometryGroup, 'photometry', AperPhotometry_Description, title='Aperture Photometry Data')
    
    for i in range(len(im_dict['startTimes'])):
        row = photometryTable.row
        row['time'] = im_dict['startTimes'][i]
        row['integration_time'] = im_dict['intTimes'][i]
        row['flux'] = fluxDict_list[i]['flux']
        row['sky'] = fluxDict_list[i]['sky']
        row['apertureRad'] = fluxDict_list[i]['apertureRad']
        row['annulusInnerRad'] = fluxDict_list[i]['annulusInnerRad']
        row['annulusOuterRad'] = fluxDict_list[i]['annulusOuterRad']
        row['flag'] = fluxDict_list[i]['flag']
        row.append()
        
    # flush the table's I/O buffer to write the data to disk
    photometryTable.flush()
    if verbose:
        print "Wrote to: "+filename
    # close the file, flush all remaining buffers
    photometryFile.close()


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
            annulus_inner - double or list of doubles of same length as self.centroid. Gives radius of inner part of sky annulus for use in "median" sky sub. Ignored if sky_sub is not "median"
            annulus_outer - double or list of doubles of same length as self.centroid. Gives radius of outer part of sky annulus for use in "median" sky sub. Ignored if sky_sub is not "median"
            interpolation - select type of interpolation to perform on image before doing photometry
                            None - to skip interpolation
                            "cubic" - to do 2D cubic interpolation
                            "linear" - to do 2D linear interpolation (default)
            Returns: Dictionary with keywords
            flux - array of flux values. [target_flux, target_sky_flux, ref0_flux, ...]
        '''

        flux = []
        sky = []
        apRadOut = []
        annInOut = []
        annOutOut = []
        flag = []

        #if sky fitting is selected for sky subtraction, do masking and produce polynomial sky image
        if sky_sub == "fit":
            print "not implemented yet! taking sky annulus instead"
            sky_sub = "median"
            #skyIm = self.fitMaskedSky(aper_radius)

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

        return {'flux': np.asarray(flux), 'sky': np.asarray(sky), 'apertureRad':np.asarray(apRadOut), 'annulusInnerRad':np.asarray(annInOut), 'annulusOuterRad':np.asarray(annOutOut), 'flag':flag}


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

    def fitMaskedSky(self, im):
        '''
        fit polynomial to sky background to estimate sky under aperture
        '''
        pass









