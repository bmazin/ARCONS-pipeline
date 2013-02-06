'''
Quick dirty code to build a tiled image of Arp 147 from the Palomar 2012 data.
JvE Jan 30 2013
'''

import numpy as np
import datetime as dt
import util.ObsFile as ObsFile
import util.utils as utils
import os.path
import hotpix.hotPixels as hotPixels
import matplotlib.pylab as mpl
import pyfits
import pickle

def makeArp147Image(obsFileList='inputFiles.list',
                    offsetLog='raster_20121209-052620-060135-edited.log',
                    wvlCalFileName='calsol_20121209-060704.h5',
                    flatCalFileName='flatsol_20121210.h5',
                    workingDir='/Users/vaneyken/Data/Palomar2012/arp147/Arp147WorkingDir/',
                    showMasks=False,
                    outputFile=''):
    
    plateScale = 0.336281230787  #0.402      #in arcsec per pix.
    imAngularWidth = 60       #Virtual image width in arcsec 
    imAngularHeight = 60      #Ditto height
    theta = 0.257             #Angle between north and pixel column direction (rad)
                          #(assume constant...!)
    slewTime = 2          #Block out this much time after each slew start (sec).
    
    wavebands = [{'color':'red', 'min':6000, 'max':8000}, #Define names and ranges (Ang) 
                 {'color':'green', 'min':5000, 'max':6000}, #for the wavebands.
                 {'color':'blue', 'min':4000, 'max':5000}, #NOTE - REDDEST BAND MUST GO FIRST.
                ]
    
    #wavebands = [{'color':'clear',   'min':4000,  'max':8000}] #Define names and ranges (Ang) 

    ##Find the index of the reddest band in 'wavebands'
    #reddestBand = np.argmax([(x['min']+x['max'])/2 for x in wavebands])
    
    raOffset = np.array([])       #To take a list of RA offsets
    decOffset = np.array([])      #Ditto for dec.
    tSlewStart = np.array([])     #For list of slew start times for each offset 
    obsStart = np.array([])       #To take list of start times for each obs. file
    obsEnd = np.array([])         #Ditto end times.
    
    #Unit vectors in RA and dec directions on virtual grid (in pixels, assume const...)
    raUnitVec = -np.array([np.cos(theta), np.sin(theta)]) / plateScale
    decUnitVec = -np.array([-np.sin(theta), np.cos(theta)]) / plateScale
    
    #Read in file name list
    #try:
    f = open(os.path.join(workingDir, obsFileList), 'r')
    obsFileNames = np.array([os.path.join(workingDir, line.strip()) for line in f.readlines()])
    #finally:
    #    print 'Problem reading file list'
    #    f.close()
    
    #Read in offset log
    try:
        f = open(os.path.join(workingDir, offsetLog), 'r')
        for eachLine in f:
            splitLine = eachLine.split()
            if len(splitLine) >= 3:             #If <3 entries in line, just skip it. (Takes care of empty lines).
                raOffset = np.append(raOffset, float(splitLine[1]))
                decOffset = np.append(decOffset, float(splitLine[2]))
                tSlewStart = np.append(tSlewStart, dt.datetime.strptime(splitLine[0],
                                                                '%Y%m%d_%H%M%S'))
    finally:
        f.close()


    #Start and end times for each offset  location on the sky
    tLocStart = tSlewStart + dt.timedelta(seconds=slewTime) #Start time of integration for each pointing
    tLocEnd = np.append(tSlewStart[1:], np.NaN)     #Last entry has no duration specified...

    #Chop off the last pointing entry from consideration since there is no way 
    #Knowing its duration
    tLocStart = tLocStart[0:-1]
    tLocEnd = tLocEnd[0:-1]
    raOffset = raOffset[0:-1]
    decOffset = decOffset[0:-1]
    nLoc = len(tLocStart)            #Number of pointings (locations) on sky

    #Get start and end times for each obs file
    for eachFileName in obsFileNames:
        obsFile = ObsFile.ObsFile(eachFileName)
        tStart = dt.datetime.strptime(obsFile.getFromHeader('utc'),
                         '%a, %d %b %Y %H:%M:%S')
        obsLength = dt.timedelta(seconds=float(obsFile.getFromHeader('exptime')))
        tEnd = tStart + obsLength
        obsStart = np.append(obsStart, tStart)
        obsEnd = np.append(obsEnd, tEnd)
        print eachFileName, tStart, tEnd

    nDetRow = obsFile.nRow  #Number of rows and columns on the detector
    nDetCol = obsFile.nCol

    print 'nDetRow, nDetCol:', nDetRow, nDetCol

    #Virtual grid width and height in pixels
    nVGridRow = np.ceil(imAngularHeight / plateScale)
    nVGridCol = np.ceil(imAngularWidth / plateScale)
    
    #Location of top left corner of a detector image if it were centered on the virtual grid.
    xCornerCen = np.floor(nVGridCol / 2 - nDetCol / 2)
    yCornerCen = np.floor(nVGridRow / 2 - nDetRow / 2)

    #Initialise 'images' dictionary to look like:
    #{'red':A, 'green':B, 'blue':C} (or whatever)
    #Where A, B, C (,...) are 3D arrays (npointings x representing an image for each pointing. 
    images = {}
    for eachWaveband in wavebands:
        images[eachWaveband['color']] = np.zeros((nLoc, nVGridRow, nVGridCol))
        images[eachWaveband['color']].fill(np.NaN)
    
    #assert 1==0
    #Loop through each raster scan location except the last (since the last has no specified duration)
    print 'Looping through pointings'
    for (iLoc, eachLocStart, eachLocEnd, eachRaOffset, eachDecOffset) in \
    zip(range(len(tLocStart)), tLocStart, tLocEnd, raOffset, decOffset):
        
        print 'Offset (RA, dec, in arcsec): ', eachRaOffset, eachDecOffset
        print 'Time: ', eachLocStart, ' to ', eachLocEnd        
        #Find which obs files have a start or end time within the time-span of the current location
        atThisLoc = ~((obsStart >= eachLocEnd) | (obsEnd <= eachLocStart))
        obsNamesAtLoc = obsFileNames[atThisLoc]     #The file names for obs files which overlap this pointing time
        obsStartAtLoc = obsStart[atThisLoc]         #The corresponding start times for those files
        obsEndAtLoc = obsEnd[atThisLoc]
                        
        #Make an image in each waveband for this location
        for eachWaveband in wavebands:
            
            print eachWaveband['color']
            #Stack all the images from all obs files who's times overlap this pointing
            im = np.zeros((nDetRow, nDetCol))    #Initialise an empty image
            totIntTime = 0                       #To take summed integration time (sec) for this pointing.
            for eachFileName, eachObsStart, eachObsEnd in zip(obsNamesAtLoc, obsStartAtLoc, obsEndAtLoc):
                print eachFileName
                obsFile = ObsFile.ObsFile(eachFileName)
                obsFile.loadFlatCalFile(os.path.join(workingDir, flatCalFileName))
                obsFile.loadWvlCalFile(os.path.join(workingDir, wvlCalFileName))
                obsFile.setWvlCutoffs(eachWaveband['min'], eachWaveband['max'])
                #Find integration start/end times in seconds relative to start of obs. file
                #Note that negative datetime objects give unexpected results - so don't go there!
                if eachLocStart < eachObsStart:
                    firstSec = 0
                else:
                    firstSec = (eachLocStart - eachObsStart).seconds      #No. of seconds from start time of obs file.
                if eachLocEnd > eachObsEnd:
                    lastSec = (eachObsEnd - eachObsStart).seconds        #Will integrate to the end of the obs. file.
                else:
                    lastSec = (eachLocEnd - eachObsStart).seconds
                integrationTime = lastSec - firstSec
                im += obsFile.getPixelCountImage(firstSec=firstSec,
                                           integrationTime=integrationTime,
                                           weighted=True)
                totIntTime+=integrationTime     #In seconds
                #assert 1==0
                del(obsFile)
            
            #Mask any bad pixels
            if eachWaveband['color'] == wavebands[0]['color']:       #If we're on the reddest band
                hpMask = hotPixels.checkInterval(image=im, fwhm=4.0)['mask']  #make a pixel mask.
                hpMask = (hpMask > 0)   #Turn into a boolean array.

            deadMask = (im == 0)                #Dead pixels
            combinedMask = hpMask | deadMask    #Always uses the hot pixel mask from the reddest (first) band.

            if showMasks:
                imCopy = np.array(im)      #Make a copy, because otherwise changes to im apparently can show up in matshow even AFTER the matshow call!!
                mpl.matshow(imCopy, vmax=np.percentile(imCopy[~np.isnan(imCopy)], 99.5))
                mpl.colorbar()
                x = np.arange(np.shape(im)[1])
                y = np.arange(np.shape(im)[0])
                xx, yy = np.meshgrid(x, y)
                mpl.scatter(xx[hpMask], yy[hpMask], c='y')
                mpl.title('#' + str(iLoc) + ', ' + str(eachLocStart) + ' to ' + str(eachLocEnd)
                         + 's, ' + str(eachWaveband['color']))
                
            im[combinedMask] = np.nan     #Set any bad pixels to NaN.
       
            #assert 1==0
            #Scale the image for the integration time:
            print 'Total integration time: '+str(totIntTime)
            im /= totIntTime
       
       
            #Determine where image should be located on the virtual grid
            gridOffset = eachRaOffset * raUnitVec + eachDecOffset * decUnitVec
            xCorner = np.round(xCornerCen + gridOffset[0])
            yCorner = np.round(yCornerCen + gridOffset[1])
            print 'X,Y corners:', xCorner, yCorner
            
            #And put it in the image stack, in the right location. All values in
            #virtual grid outside the current image will remain NaN for this pointing.
            images[eachWaveband['color']][iLoc, yCorner:yCorner + nDetRow, xCorner:xCorner + nDetCol] = im
            
            #assert 1==0
    #Now have a stack of images in 'images'

    finalImage = {}     #Will end up as dictionary with one image for each color.

    print 'Combining images'    
    for eachWaveband in wavebands:
        
        print eachWaveband['color']
        #Make numpy masked array of the stack where NaNs are all masked
        masked = np.ma.masked_invalid(images[eachWaveband['color']])
    
        #Median combine the images (leaving out masked pixels)
        finalImage[eachWaveband['color']] = np.ma.median(masked, axis=0).data
        #assert 1==0
    print 'Done....'
    #assert 1==0
    
    result = {'final':finalImage, 'stack':images}
    
    if type(outputFile) is str and len(outputFile) > 0:
        try:
            f = open(os.path.join(workingDir, outputFile), 'wb')
            pickle.dump(result, f)
        except OSError:
            print 'Could not write results to output file.' 
        finally:
            f.close()
        
    return result

    



def makeRGBimage(result, fitsBaseName='Arp147-',
                  workingDir='/Users/vaneyken/Data/Palomar2012/arp147/Arp147WorkingDir/',
                  interpolate=False, medianSmooth=False):
    '''
    Take a 'result' (as returned by makeArp147Image()) and make an RGB image out
    of it. Also create fits files of the three components.

    INPUTS:
        result: as returned by makeArp147Image()
        fitsBaseName: fits files are output as fistBaseName+[color]+'.fits'
    
    '''
    
    rgb = result['final']       #Convenient shorthand
    
    #Write to fits files
    for color in rgb.keys():
        fileName = fitsBaseName + color + '.fits'
        component = rgb[color]
        if interpolate:
            component[component == 0] = np.nan
            component = utils.replaceNaN(component, iterate=False)
        if medianSmooth:
            component = utils.median_filterNaN(component, size=3)
            
        hdu = pyfits.PrimaryHDU(component)
        hdu.writeto(os.path.join(workingDir, fileName))
        





if __name__ == "__main__":
    
    out1 = makeArp147Image(showMasks=True,
                               offsetLog='raster_20121209-052620-060135-edited-short.log')
    