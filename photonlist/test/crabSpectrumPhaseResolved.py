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

def func(wvls, K, alph):  #this is the power law fit
    return K*(wvls/6000)**-(alph+2)

def Image():
   #output = open('/home/kingrice/ARCONS-pipeline/examples/crabN/stackedImageCrab.pkl', 'r')
   output = open('/home/kingrice/ARCONS-pipeline/examples/crabN/CrabStack_phaseWvl.pkl', 'rb')
   results = pickle.load(output)
   vim = results['vim']  #pulls out the dictionary containing all of the relevant information for virtual image object
   #this is the info from the dictionary:
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

   print 'ob.image', np.shape(ob.image)

   
   #print 'ob.image', ob.image, np.shape(ob.image)
   #print 'ob.expTimeWeights', ob.expTimeWeights, np.shape(ob.expTimeWeights)
   #image = np.sum(ob.image*ob.expTimeWeights, axis = 2)
   
   image = ob.image/ob.effIntTimes #divides image cube by effective integration times to get a veiwable image. This will be 4D here.
   nanspot = np.isnan(image) #sets the NANS of the image = 0 for the display
   image[nanspot] = 0
   ob.display(DDArr = True) #sums across the wavlength and phase dimension to get a 2D image for the display.
   mpl.show()
   
   #spectrum, spectrumP, spectrumIP, phaseProfile, wvlBinEdges, phaseBinEdges, innerAp = ob.getApertureSpectrum(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.0005, radius2 = .0008, degrees=True, error = False, phase = True, crab = True) #this calls the spectrum code
   print 'getting phase of crab now' 
   spectrum, spectrumP, spectrumIP, phaseProfile, wvlBinEdges, phaseBinEdges, innerAp = ob.getAperturePhaseCrab(cenRA=83.6330495, cenDec=22.0144933, nPixRA=250, nPixDec=250, radius1=.0005, degrees=True, phase = True, crab = True)
    
   #summed_arrayIn, wvlBinEdges, errIn, innerAp = ob.getApertureSpectrum(cenRA=83.6335573, cenDec=22.0142943, nPixRA=250, nPixDec=250, radius1=.00015, radius2 = .0002, degrees=True, error = True, offSet = True, radRA = 83.6296510, radDec = 22.0145475, radius3 = .0005) #object spectrum
   
   print 'wavelengths:', wvlBinEdges[16:25] #this is the desired wavlength range for the crab spectrum it may be worth checking.

   #test
   #summed_arrayIn, wvlBinEdges, errIn, innerAp = ob.getApertureSpectrum(cenRA=83.6308, cenDec=22.0151, nPixRA=250, nPixDec=250, radius1=.0007, radius2 = .0014, degrees=True, error = True) 
   #errIn = errIn/np.mean(summed_arrayIn)
   #summed_arrayIn = summed_arrayIn/np.mean(summed_arrayIn)
   
   c=3.00E18 #Angs/s
   h=6.626E-27 #erg*s
   
   #finds the center of the wavelength bins     
   nWvlBins = len(wvlBinEdges) - 1  
         
   wvls = np.empty((nWvlBins), dtype = float)
   for n in xrange(nWvlBins):
      binsize = wvlBinEdges[n+1] - wvlBinEdges[n]
      wvls[n] = (wvlBinEdges[n] + (binsize/2.0))
   
   #finds the center of the phase bins
   nphaseBins = len(phaseBinEdges) - 1
   
   phases = np.empty((nphaseBins), dtype = float)
   for n in xrange(nphaseBins):
      pbinsize = phaseBinEdges[n+1]-phaseBinEdges[n]
      phases[n] = (phaseBinEdges[n] + (pbinsize/2.0))


   #this is where the dereddening happens
   Alam = (-7.51*np.log10(wvls/10000) + 1.15)*0.51 #find the extinction magnitude for each wavelength bin in terms of magnitude     
   
   AlamF = spectrum*10**(-Alam/2.5)  #converts to flux for entire spectrum of all time
   AlamFP = spectrumP*10**(-Alam/2.5) #this is specifically for the pulse spectrum
   AlamFIP = spectrumIP*10**(-Alam/2.5) #interpulse spectrum
   #AlamFErr = errIn*10**(-Alam/2.5)   

   spectrum = spectrum - AlamF  #subtracts from the flux values we observed
   spectrumP = spectrumP - AlamFP #main pulse
   spectrumIP = spectrumIP - AlamFIP #interPulse
  
   #errIn = np.sqrt(errIn**2+ AlamFErr**2)
   
   
   diam = 510.55 #5 meter telescope in centimeters.
   area = np.pi * ((diam/2.0)**2 -(183/2.0)**2) #secondary obstruction diameter 1.83m
         
   spectrum = spectrum*(h*c/wvls)*(1/area)*(1./float(ob.totExpTime)) #(count/Angs)*(erg*sec/count)*(Angs/s)/(Angs)/(cm^2)/(s) = erg/s/Angs/cm^2

   spectrumP = spectrumP*(h*c/wvls)*(1./area)*(50./6.)*(1./float(ob.totExpTime)) #same units as above only now scaled by number of phase bins used
   spectrumIP = spectrumIP*(h*c/wvls)*(1./area)*(50./6.)*(1./float(ob.totExpTime))
   #errIn = errIn*(h*c/wvls)
   #errIn = errIn*(h*c/wvls)*(1/area)   
   #print 'errIn', errIn, np.shape(errIn)

   print 'spectrum', spectrum, np.shape(spectrum)
   print 'spectrumP', spectrumP, np.shape(spectrumP)
   print 'spectrumIP', spectrumIP, np.shape(spectrumIP)
   #note that in this particular case we are plotting the spectrum on a log scale. In order to convert te errors do the following.

   #errlog = (1/(summed_arrayIn*np.log(10)))*errIn

   '''
   K = 5.9*10**-15
   lam0 = 6000
   alph = .2
   
   #func = K*(wvls/lam0)**-(alph+2)
   params = [K,lam0,alph]
   errs = errIn
   quiet = True
   parinfo = [{'n':0,'value':params[0],'limits':[1.9*10**-15,10*10**-15],'limited':[True,True],'fixed':False,'parname':"k",'error':0},{'n':1,'value':params[1],'limits':[3000,9000],'limited':[True,True],'fixed':False,'parname':"lam0",'error':0}, {'n':2,'value':params[2],'limits':[-.5,.5],'limited':[True,True],'fixed':False,'parname':"alph",'error':0}]
   
   fa = {'x':wvls,'y':summed_arrayIn,'err':errIn}

   m = mpfit.mpfit(func, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
   
   if m.status<=0:
      print m.status, m.errmsg

   mpp = m.params
   mpperr = m.perror
   

   for k,p in enumerate(mpp):
      parinfo[k]['value'] = p
      parinfo[k]['error'] = mpperr[k]
      if k==0: K = p
      if k==1: lam0 = p
      if k==2: alph = p

   fit = K*(wvls/lam0)**-(alph+2)
   '''
   
   print 'wvls', wvls, np.shape(wvls)
   x0 = np.array([5.9*10**-15,.2]) #initial guesses for fit paramaters for main pulse
   x1 = np.array([1.9*10**-15,.2]) #initial guesses for fit paramaters for main pulse
   

   #this makes the fit for the various spectra
   popt, pcov =  optimization.curve_fit(func, wvls[15:28], spectrum[15:28], x0) #all
   poptP, pcovP = optimization.curve_fit(func, wvls[15:28], spectrumP[15:28], x0) # main pulse
   poptIP, pcovIP = optimization.curve_fit(func, wvls[15:28], spectrumIP[15:28], x1) #interpulse
   
   print popt
   print poptP
   print poptIP
    
   #calculates the fit bases on the fit parameters
   fit = popt[0]*(wvls[15:28]/6000)**-(popt[1]+2) 
   fitP = poptP[0]*(wvls[15:28]/6000)**-(poptP[1]+2)  # main pulse
   fitIP = poptIP[0]*(wvls[15:28]/6000)**-(poptIP[1]+2) #interpulse
   
   '''
   plotArray(np.sum(image[:,:,0:5],axis=2))
   plotArray(np.sum(image[:,:,5:10],axis=2))
   plotArray(np.sum(image[:,:,10:15],axis=2))
   plotArray(np.sum(image[:,:,15:20],axis=2))
   plotArray(np.sum(image[:,:,20:25],axis=2))
   plotArray(np.sum(image[:,:,25:31],axis=2))
   '''
   #plotArray(image)
   
   ob.display(DDArr = True)
   mpl.show()
   
   mpl.figure(1) #for all phases
   mpl.step(wvls,np.log10(spectrum), color = 'r', where = 'pre')
   #mpl.errorbar(wvls,np.log10(summed_arrayIn), yerr = errlog, ecolor = 'b', fmt = None)
   mpl.plot(wvls[15:28], np.log10(fit), color = 'b')
   mpl.xlim(wvls[15], 9000)
       
   mpl.figure(2) #this plot is the main pulse and interpulse
   mpl.step(wvls,np.log10(spectrumP), color = 'r', where = 'pre')
   #mpl.errorbar(wvls,np.log10(summed_arrayIn), yerr = errlog, ecolor = 'b', fmt = None)
   mpl.plot(wvls[15:28], np.log10(fitP), color = 'b')
   mpl.xlim(wvls[15], 9000)
    
   mpl.step(wvls,np.log10(spectrumIP), color = 'g', where = 'pre')
   #mpl.errorbar(wvls,np.log10(summed_arrayIn), yerr = errlog, ecolor = 'b', fmt = None)
   mpl.plot(wvls[15:28], np.log10(fitIP), color = 'b')
   mpl.xlim(wvls[15], 9000)     

   mpl.figure(3) #this is the phase profile
   mpl.step(phases, phaseProfile, color = 'b', where = 'pre')
   '''
   mpl.plot(wvMid, np.log10(FlamMid), color = 'b')
   mpl.plot(wvHigh, np.log10(FlamHigh), color = 'b')
   mpl.plot(wvLow, np.log10(FlamLow), color = 'b')
   mpl.plot(wvPrev, np.log10(FlamPrev), color = 'g')
   #mpl.errorbar(wvlBinEdges[1:],np.log(summed_arrayIn), yerr = errIn, ecolor = 'b', fmt = None)
   mpl.xlim(4000, 9000)
   mpl.ylim(-14.6, -13.9) 
   
   mpl.figure(2)
   '''
   mpl.show()

if __name__=='__main__':
   Image()
