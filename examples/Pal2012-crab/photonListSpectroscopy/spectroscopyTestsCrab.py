import numpy as np
import warnings
import pickle
import os.path
import glob
import scipy.stats
import numpy as np
#from astropy import coordinates as coord
import matplotlib.pyplot as mpl
import photonlist.photlist as pl
from photonlist.RADecImage import RADecImage
from util.FileName import FileName
from util import utils
from util.popup import PopUp, plotArray
import scipy.optimize as optimization
from util import mpfit

"""
written by Neil King Rice - October, 2014

This code is designed to test the three aperture methods in the RADecImage Class written for analysis on data cubes generated from photon lists. the data used here is for the Crab Pulsar 
The three methods are:
    getApertureSpectrum
    getAperturePhase
    getPhaseSpectrum

This code also includes some fitting routines.

"""

def func(wvls, K, alph):  #this is the power law fit
    return K*(wvls/6000)**-(alph+2)

def backgroundAnalysis():
    output = open('/home/kingrice/ARCONS-pipeline/examples/crabN/CrabStack_phaseWvl.pkl', 'rb')
    results = pickle.load(output)
    vim = results['vim']  #pulls out the dictionary containing all of the relevant information for virtual image object
   
    #this is the info included in the dictionary:
    ob = RADecImage()
    ob.phaseEdges = vim['phase']   
    ob.wvlEdges = vim['wvl']
    ob.gridDec = vim['Dec']
    ob.gridRA = vim['RA']
    ob.image = vim['Im']
    ob.effIntTimes = vim['effIntT']
    ob.totExpTime = vim['totExpT']
    ob.expTimeWeights = vim['timeWeights']
    ob.vExpTimesStack = vim['vTimeStacks']
    ob.imageIsLoaded = vim['loaded']
    output.close()

    ob.display(DDArr = True) #sums across the wavlength and phase dimension to get a 2D image for the display.
    mpl.show()

    '''
    Here is the fir4st step in the analysis. We want to obtain a sky subtracted spectrum for lets say 5 different scenarios where
    the same aperture is used but six different anuli are used in addition to the original annulus. This should be no problem.
    '''

    spectrumOG, skySubtractionOG, innerApOG, skyMaskOG, wvlBinEdgesOG, apAreaOG, anAreaOG = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, radius2 = .000770, degrees=True) #OG

    spectrumT1, skySubtractionT1, innerApT1, skyMaskT1, wvlBinEdgesT1, apAreaT1, anAreaT1 = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, radius2 = .000530, degrees=True, offSet=True, radRA=83.6301176, radDec=22.0143848, radius3=.000770) #T1
    
    spectrumT2, skySubtractionT2, innerApT2, skyMaskT2, wvlBinEdgesT2, apAreaT2, anAreaT2 = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, radius2 = .000530, degrees=True, offSet=True, radRA=83.6317452, radDec=22.0112381, radius3=.000770) #T2

    spectrumT3, skySubtractionT3, innerApT3, skyMaskT3, wvlBinEdgesT3, apAreaT3, anAreaT3 = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, radius2 = .000530, degrees=True, offSet=True, radRA=83.6338792, radDec=22.0119614, radius3=.000770) #T3

    spectrumT4, skySubtractionT4, innerApT4, skyMaskT4, wvlBinEdgesT4, apAreaT4, anAreaT4 = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, radius2 = .000530, degrees=True, offSet=True, radRA=83.6350004, radDec=22.0135529, radius3=.000770) #T4

    spectrumT5, skySubtractionT5, innerApT5, skyMaskT5, wvlBinEdgesT5, apAreaT5, anAreaT5 = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, radius2 = .000530, degrees=True, offSet=True, radRA=83.6341686, radDec=22.0179293, radius3=.000770) #T5

    spectrumT6, skySubtractionT6, innerApT6, skyMaskT6, wvlBinEdgesT6, apAreaT6, anAreaT6 = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, radius2 = .000530, degrees=True, offSet=True, radRA=83.6350004, radDec=22.0161571, radius3=.000770) #T6

    
    print 'apArea = ', apAreaOG, ' (arcseconds^2)'
    print 'anArea = ', anAreaOG, ' (arcseconds^2)'    
    
    print '---determining wavelengthBin centers---'

    nWvlBins = len(wvlBinEdgesOG) - 1  
         
    wvls = np.empty((nWvlBins), dtype = float)
    for n in xrange(nWvlBins):
        binsize = wvlBinEdgesOG[n+1] - wvlBinEdgesOG[n]
        wvls[n] = (wvlBinEdgesOG[n] + (binsize/2.0))

    print '---gathering conversion factors---'
    
    c=2.998E18 #Angs/s
    h=6.626E-27 #erg*s
    diam = 510.55 #5 meter telescope in centimeters.
    area = np.pi * ((diam/2.0)**2 -(183/2.0)**2) #secondary obstruction diameter 1.83m 

    print '---dereddening spectra---'
   
    #this is where the dereddening happens
    Alam = (-7.51*np.log10(wvls/10000) + 1.15)*0.52 #find the extinction magnitude for each wavelength bin in terms of magnitude

    spectrumOG = spectrumOG*10**(Alam/2.5)
    spectrumT1 = spectrumT1*10**(Alam/2.5)
    spectrumT2 = spectrumT2*10**(Alam/2.5) 
    spectrumT3 = spectrumT3*10**(Alam/2.5) 
    spectrumT4 = spectrumT4*10**(Alam/2.5) 
    spectrumT5 = spectrumT5*10**(Alam/2.5) 
    spectrumT6 = spectrumT6*10**(Alam/2.5)

    skySubtractionOG = skySubtractionOG*10**(Alam/2.5)
    skySubtractionT1 = skySubtractionT1*10**(Alam/2.5)
    skySubtractionT2 = skySubtractionT2*10**(Alam/2.5)
    skySubtractionT3 = skySubtractionT3*10**(Alam/2.5)
    skySubtractionT4 = skySubtractionT4*10**(Alam/2.5)
    skySubtractionT5 = skySubtractionT5*10**(Alam/2.5)
    skySubtractionT6 = skySubtractionT6*10**(Alam/2.5)   

    print '---converting flux units to erg/s/cm^2/A---'

    spectrumOG = spectrumOG*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    spectrumT1 = spectrumT1*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    spectrumT2 = spectrumT2*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    spectrumT3 = spectrumT3*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    spectrumT4 = spectrumT4*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    spectrumT5 = spectrumT5*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    spectrumT6 = spectrumT6*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    
    skySubtractionOG = skySubtractionOG*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    skySubtractionT1 = skySubtractionT1*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    skySubtractionT2 = skySubtractionT2*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    skySubtractionT3 = skySubtractionT3*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    skySubtractionT4 = skySubtractionT4*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    skySubtractionT5 = skySubtractionT5*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    skySubtractionT6 = skySubtractionT6*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    

    print '---plotting photometry points for comparison---'
   
    #photPoints = np.array([17.23,16.66,16.17,15.65]) 
    photPoints = np.array([17.22,16.64,16.14,15.61])
    photPointsErr = np.array([.02, .02, .01, .01])
    photPointsFlux = np.array([4260.,3640.,3080.,2550.]) #these are in janskys
    centerWvl = np.array([4400.,5500.,6400.,7900.])   #these are in Angstrom

    AlamPP = (-7.51*np.log10(centerWvl/10000.) + 1.15)*0.51

    photPointsFlux = photPointsFlux*(10**(AlamPP/2.5)) #dereddened
    photPointsErr = photPointsErr*(10**(AlamPP/2.5))
    # 1Jy = 10^-23 erg/sec/cm^2/Hz   
    
    photPointsFlux = photPointsFlux*(10**-23)  #now in units of erg/sec/cm^2/Hz
    photPointsErr = photPointsErr*(10**-23)

    photPointsFlux = photPointsFlux*(c/(centerWvl**2)) #now in units of erg/s/cm^2/A
    photPointsErr = photPointsErr*(c/(centerWvl**2))

    ppfFinal = photPointsFlux*(10**(-photPoints/2.5))  #now same units as plot
    ppEFinal = photPointsErr*(10**(-photPoints/2.5))


    image = (ob.image*ob.expTimeWeights) 
    image = np.sum(image, axis = 3)
    image = np.sum(image, axis = 2)
    image[np.where(skyMaskOG==0)] = 0
    image[np.where(skyMaskT1==0)] = 0
    image[np.where(skyMaskT2==0)] = 0
    image[np.where(skyMaskT3==0)] = 0
    image[np.where(skyMaskT4==0)] = 0
    image[np.where(skyMaskT5==0)] = 0
    image[np.where(skyMaskT6==0)] = 0
    #image[np.where(innerApOG==0)] = 0
    #image[np.where(innerApOG==0)] = 0
    #image[np.where(innerApOG==0)] = 0    
    
    ob.display(image = image) #sums across the wavlength and phase dimension to get a 2D image for the display.
    mpl.show()

    mpl.figure(1)
    mpl.step(wvls, spectrumOG*(10**14), color = 'r', where = 'mid')
    mpl.step(wvls, spectrumT1*(10**14), color = 'b', where = 'mid')
    mpl.step(wvls, spectrumT2*(10**14), color = 'g', where = 'mid')
    mpl.step(wvls, spectrumT3*(10**14), color = 'y', where = 'mid')
    mpl.step(wvls, spectrumT4*(10**14), color = 'c', where = 'mid')
    mpl.step(wvls, spectrumT5*(10**14), color = 'k', where = 'mid')
    mpl.step(wvls, spectrumT6*(10**14), color = 'm', where = 'mid')
    mpl.plot(centerWvl, ppfFinal*(10**14), 'bo')
    mpl.xlim(4000,10000)
    mpl.ylim(0,1)

    mpl.figure(2)
    mpl.step(wvls, skySubtractionOG*(10**16), color = 'r', where = 'mid')
    mpl.step(wvls, skySubtractionT1*(10**16), color = 'b', where = 'mid')
    mpl.step(wvls, skySubtractionT2*(10**16), color = 'g', where = 'mid')
    mpl.step(wvls, skySubtractionT3*(10**16), color = 'y', where = 'mid')
    mpl.step(wvls, skySubtractionT4*(10**16), color = 'c', where = 'mid')
    mpl.step(wvls, skySubtractionT5*(10**16), color = 'k', where = 'mid')
    mpl.step(wvls, skySubtractionT6*(10**16), color = 'm', where = 'mid')
    mpl.xlim(4000,10000)
    mpl.ylim(0,1)


    mpl.show()
    ###############################################################################################################
    '''
    This next step is to demonstrate the functionality of the getApertureSpectrum when phase is included
    '''


    spectrumAlt, phaseSpectrumAlt, skySubtractionAlt, innerAppAlt, skyMaskAlt, wvlBinEdgesAlt, phaseBinEdgesAlt, apAreaAlt, anAreaAlt = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, radius2 = .000770, degrees=True, phase = True)
        
    skySubtractionAlt = np.sum(skySubtractionAlt, axis = 1)
    phaseProFileAlt = np.sum(phaseSpectrumAlt, axis = 0)


    print '---determining wavelengthBin centers---'

    nWvlBins = len(wvlBinEdgesAlt) - 1  
         
    wvls = np.empty((nWvlBins), dtype = float)
    for n in xrange(nWvlBins):
        binsize = wvlBinEdgesAlt[n+1] - wvlBinEdgesAlt[n]
        wvls[n] = (wvlBinEdgesAlt[n] + (binsize/2.0))

    print '---determining phase bin centers---'

    nphaseBins = len(phaseBinEdgesAlt) - 1
   
    phases = np.empty((nphaseBins), dtype = float)
    for n in xrange(nphaseBins):
        pbinsize = phaseBinEdgesAlt[n+1]-phaseBinEdgesAlt[n]
        phases[n] = (phaseBinEdgesAlt[n] + (pbinsize/2.0))

    print '---gathering conversion factors---'
    
    c=2.998E18 #Angs/s
    h=6.626E-27 #erg*s
    diam = 510.55 #5 meter telescope in centimeters.
    area = np.pi * ((diam/2.0)**2 -(183/2.0)**2) #secondary obstruction diameter 1.83m 

    print '---dereddening spectra---'
   
    #this is where the dereddening happens
    Alam = (-7.51*np.log10(wvls/10000) + 1.15)*0.52 #find the extinction magnitude for each wavelength bin in terms of magnitude

    spectrumAlt = spectrumAlt*10**(Alam/2.5)
    skySubtractionAlt = skySubtractionAlt*10**(Alam/2.5)

    print '---converting flux units to erg/s/cm^2/A---'

    spectrumAlt = spectrumAlt*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    skySubtractionAlt = skySubtractionAlt*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))
    

    mpl.figure(3)
    mpl.step(wvls, spectrumAlt*(10**14), color = 'r', where = 'mid')
    mpl.plot(centerWvl, ppfFinal*(10**14), 'bo')
    mpl.xlim(4000,10000)
    mpl.ylim(0,1)
    
    mpl.figure(4)
    mpl.step(wvls, skySubtractionAlt*(10**16), color = 'r', where = 'mid')
    mpl.xlim(4000,10000)
    mpl.ylim(0,1)    

    mpl.figure(5)
    mpl.step(phases, phaseProFileAlt, color = 'k', where = 'mid')

    mpl.show()

    #############################################################################################################################
    '''
    the previous step was to get an idea of how the sky background changes from place to place in the sky. Now we want to compare the
    off pulse flux to the various sky annuli. In adittion we should now like to plot the off pulse region 
    '''
    
    phaseWvlCount, wvlBinEdges, phaseBinEdges, innerAp, apArea = ob.getAperturePhase(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.000530, degrees=True, spectrum = True)

    print np.shape(phaseWvlCount)

    phaseProfile = np.sum(phaseWvlCount, axis = 1)
    phaseProfile = np.sum(phaseProfile, axis = 0)
    
    
    print np.shape(phaseProfile)

    
    print '---determining phase bin centers---'

    nphaseBins = len(phaseBinEdges) - 1
   
    phases = np.empty((nphaseBins), dtype = float)
    for n in xrange(nphaseBins):
        pbinsize = phaseBinEdges[n+1]-phaseBinEdges[n]
        phases[n] = (phaseBinEdges[n] + (pbinsize/2.0))
    
    print np.shape(phases)    

    print '---grouping spectrum by phases---'
    
    #spectrumAll = ob.getPhaseSpectrum(countArray=phaseWvlCount, lowLim=0., highLim=1.)
    
    #spectrumAll = np.sum(phaseWvlCount, axis = 2)
    #spectrumAll = np.sum(spectrumAll, axis = 0)
    #for i in range(len(spectrumAll)):
    #    spectrumAll[i]/=(wvlBinEdges[i+1]-wvlBinEdges[i])


    spectrumP = ob.getPhaseSpectrum(countArray=phaseWvlCount, lowLim=.6, highLim=.74)  

    spectrumIP = ob.getPhaseSpectrum(countArray=phaseWvlCount, lowLim=.02, highLim=.16) 

    spectrumNeb1 = ob.getPhaseSpectrum(countArray=phaseWvlCount, lowLim=.4, highLim=.54, background = True, med = True) 

    spectrumNeb2 = ob.getPhaseSpectrum(countArray=phaseWvlCount, lowLim=.4, highLim=.54, background = True)

    print '---gathering conversion factors---'
    
    spectrumAll = np.sum(phaseWvlCount, axis = 2)

    for i in range(len(spectrumAll)):
        for j in range(len(spectrumAll[i])):

            spectrumAll[i][j]/=(wvlBinEdges[j+1]-wvlBinEdges[j])

    for i in range(len(spectrumNeb1)):
        spectrumNeb1[i]/=(wvlBinEdges[i+1]-wvlBinEdges[i])
    
    for i in range(len(spectrumP)):
        for j in range(len(spectrumP[i])):
            spectrumP[i][j]/=(wvlBinEdges[j+1]-wvlBinEdges[j])
 
    for i in range(len(spectrumIP)):
        for j in range(len(spectrumIP[i])):
            spectrumIP[i][j]/=(wvlBinEdges[j+1]-wvlBinEdges[j])

    for i in range(len(spectrumNeb2)):
        for j in range(len(spectrumNeb2[i])):
            spectrumNeb2[i][j]/=(wvlBinEdges[j+1]-wvlBinEdges[j])


    c=3.00E18 #Angs/s
    h=6.626E-27 #erg*s
    diam = 510.55 #5 meter telescope in centimeters.
    area = np.pi * ((diam/2.0)**2 -(183/2.0)**2) #secondary obstruction diameter 1.83m 

    print '---dereddening spectra---'
   
    #this is where the dereddening happens
    Alam = (-7.51*np.log10(wvls/10000) + 1.15)*0.51 #find the extinction magnitude for each wavelength bin in terms of magnitude     
   
    
    spectrumAll = spectrumAll*10**(Alam/2.5)
    
    

    spectrumP = spectrumP*10**(Alam/2.5) #this is specifically for the pulse spectrum
    spectrumIP = spectrumIP*10**(Alam/2.5) #interpulse spectrum
    spectrumNeb1 = spectrumNeb1*10**(Alam/2.5) #deadtime betweeen IP and P
    spectrumNeb2 = spectrumNeb2*10**(Alam/2.5) #deadtime betweeen P and IP 
   
    print '---converting flux units to erg/s/cm^2/A---'
    
    spectrumAll = spectrumAll*(h*c/wvls)*(1./area)*(1./float(ob.totExpTime))

    spectrumP = spectrumP*(h*c/wvls)*(1./area)*(50./7.)*(1./float(ob.totExpTime))
    spectrumIP = spectrumIP*(h*c/wvls)*(1./area)*(50./7.)*(1./float(ob.totExpTime))
    spectrumNeb1 = spectrumNeb1*(h*c/wvls)*(1./area)*(50./7.)*(1./float(ob.totExpTime))
    spectrumNeb2 = spectrumNeb2*(h*c/wvls)*(1./area)*(50./7.)*(1./float(ob.totExpTime))
    
    for i in range(len(spectrumAll)):
        spectrumAll[i] = spectrumAll[i] - spectrumNeb1
    
    #spectrumAllCorr = spectrumAll - spectrumNeb1 
    #spectrumPCorr = spectrumP - spectrumNeb1
    #spectrumIPCorr = spectrumIP - spectrumNeb2

    for i in range(len(spectrumP)):
        spectrumP[i] = spectrumP[i] - spectrumNeb2[i]
        
    for i in range(len(spectrumIP)):
        spectrumIP[i] = spectrumIP[i] -spectrumNeb2[i]

    spectrumAll = np.sum(spectrumAll, axis = 0)
    spectrumIP = np.sum(spectrumIP, axis = 0)
    spectrumP = np.sum(spectrumP, axis = 0)
    
    mpl.figure(6)
    mpl.step(phases, phaseProfile, color = 'k', where = 'mid')
    
    mpl.figure(7)
    mpl.step(wvls, spectrumAll*(10**14), color = 'r', where = 'mid')
    mpl.step(wvls, spectrumP*(10**14), color = 'b', where = 'mid')
    mpl.step(wvls, spectrumIP*(10**14), color = 'g', where = 'mid')
    mpl.plot(centerWvl, ppfFinal*(10**14), 'bo')
    mpl.errorbar(centerWvl, ppfFinal*(10**14), yerr = ppEFinal*(10**14), fmt = 'k.') 
    mpl.xlim(4000,10000)
    mpl.ylim(0,2.5)

    mpl.figure(8)
    mpl.step(wvls, spectrumNeb1*(10**16), color = 'r', where = 'mid')
    #mpl.step(wvls, spectrumNeb2*(10**16), color = 'b', where = 'mid')
    mpl.xlim(4000,10000)
    mpl.ylim(0,1)

    mpl.figure(9)
    mpl.step(wvls, spectrumNeb2[0]*(10**16), color = 'b', where = 'mid')
    mpl.step(wvls, spectrumNeb2[1]*(10**16), color = 'b', where = 'mid')
    mpl.step(wvls, spectrumNeb2[2]*(10**16), color = 'b', where = 'mid')
    mpl.step(wvls, spectrumNeb2[3]*(10**16), color = 'b', where = 'mid')
    mpl.step(wvls, spectrumNeb2[50]*(10**16), color = 'b', where = 'mid')
    mpl.step(wvls, spectrumNeb2[100]*(10**16), color = 'b', where = 'mid')
    mpl.step(wvls, spectrumNeb2[-1]*(10**16), color = 'b', where = 'mid')
    mpl.xlim(4000,10000)
    mpl.ylim(0,1)
    

    mpl.show()
    

    '''
    Now we apply our Fits
    '''
    spectrumNeb2 = np.sum(spectrumNeb2, axis = 0)

    print '---determining initial gusses for fit parameters---'
    x0 = np.array([5.9*10**-15,.2]) #initial guesses for fit paramaters for main pulse
    x1 = np.array([1.9*10**-15,.2]) #initial guesses for fit paramaters for interPulse
    x2 = np.array([3.4*10**-15,-.4]) #initial guesses for fit paramaters for nebula
   
    #this makes the fit for the various spectra
    print '---determining fits---'
   
    #poptP, pcovP = optimization.curve_fit(func, wvls[15:28], spectrumP[15:28], x0) # main pulse
    #poptIP, pcovIP = optimization.curve_fit(func, wvls[15:28], spectrumIP[15:28], x1) #interpulse
    #poptNeb1, pcovNeb1 = optimization.curve_fit(func, wvls[15:28], (spectrumNeb1[15:28]/apArea), x2) # main pulse
    poptNeb2, pcovNeb2 = optimization.curve_fit(func, wvls[15:28], (spectrumNeb2[15:28]/apArea), x2) #interpulse
    #poptNebMed, pcovNebMed = optimization.curve_fit(func, wvls, (nebMedian[15:28]/apArea), x2) #interpulse
    poptP, pcovP = optimization.curve_fit(func, wvls[15:28], spectrumP[15:28], x0) # main pulse
    poptIP, pcovIP = optimization.curve_fit(func, wvls[15:28], spectrumIP[15:28], x1) #interpulse

    print 'fit prameters: '
    #print poptP
    #print poptIP
    #print 'Nebula1/apArea - ', poptNeb1
    print 'Nebula2/apArea - ', poptNeb2
    #print 'Median Nebula/apArea', poptNebMed
    print 'Main Pulse - ', poptP
    print 'Interpulse - ', poptIP
    
    #fitP = poptP[0]*(wvls[15:28]/6000)**-(poptP[1]+2)  # main pulse
    #fitIP = poptIP[0]*(wvls[15:28]/6000)**-(poptIP[1]+2) #interpulse
    #fitNeb1 = poptNeb1[0]*(wvls[15:28]/6000)**-(poptNeb1[1]+2)  # nebula
    fitNeb2 = poptNeb2[0]*(wvls[15:28]/6000)**-(poptNeb2[1]+2) #nebula
    #fitNebMed = poptNebMed[0]*(wvls/6000)**-(poptNebMed[1]+2) #nebula
    fitP = poptP[0]*(wvls[15:28]/6000)**-(poptP[1]+2)  # main pulse
    fitIP = poptIP[0]*(wvls[15:28]/6000)**-(poptIP[1]+2) #interpulse

    print 'wvls[15:28]', wvls[15:28]
    print 'wvls[12:28]', wvls[12:28]
    
    mpl.figure(10)
    mpl.step(wvls, np.log10(spectrumP), color = 'b', where = 'mid')
    mpl.step(wvls, np.log10(spectrumIP), color = 'g', where = 'mid')
    mpl.plot(wvls[15:28], np.log10(fitP), color = 'k')
    mpl.plot(wvls[15:28], np.log10(fitIP), color = 'k')
    
    mpl.xlim(4000,10000)

    mpl.figure(11)
    mpl.step(wvls, np.log10(spectrumNeb2/apArea), color = 'b', where = 'mid')
    mpl.plot(wvls[15:28], np.log10(fitNeb2), color = 'k')
    
    mpl.xlim(4000,10000)
    
    mpl.show()

    
if __name__=='__main__':
    backgroundAnalysis()   
