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
from timing.photonTiming import *

def writeToPHOEBEFile(fname, objName, passband, times, mags, errors):
    '''
    function to output data to PHOEBE formatted txt file
    '''
    f = open(fname,'w')
    f.write('# Object:\t%s\n'%objName)
    f.write('# Passband:\t%s\n'%passband)
    f.write('# Source:\tARCONS, Palomar Observatory\n')
    f.write('# BJD\t%smag\tError\n'%passband)
    for i in xrange(len(times)):
        f.write('%.6f\t%.3f\t%.3f\n'%(times[i],mags[i],errors[i]))
    f.close()
    print "Finished writing to PHOEBE format file"
    return


if len(sys.argv)>1:
    filt = sys.argv[1]
else:
    filt = None

objName = '1SWASP_J0002'
path = '/Scratch/DisplayStack/PAL2014/1SWASP_J0002' #J0002 for my object, J2210 for Alex's object
verbose=True
showPlot=True

intTime=30
startWl = 3000
endWl = 13000
hp=True
flat = True
flux = True

#outputFN = '%ss_%s-%s'%(intTime,startWl, endWl)
outputFN = '%ss_%s'%(intTime,filt)

if flat==True:
    outputFN+='_flat'
if flux==True:
    outputFN+='_flux'
if hp==True:
    outputFN+='_hp'
#outputFN+='.h5'

identifier=outputFN
print outputFN


LC=LightCurve(fileID=identifier,path=path,targetName=None,run=None,verbose=True,showPlot=False)

'''
LC.loadImageStack()
LC.loadAllCentroidFiles()
print LC.centroids
print LC.flags
#print "Finished loading centroid files"
LC.makeLightCurve(photometryType='aper', interpolation='cubic')


'''


LC.loadLightCurve(photometryFilename='',photometryType='aper')
photometryDict=LC.photometry_dict

#get timestamp of bin centers in UNIX time
time = photometryDict['startTimes']
expTime = photometryDict['intTimes']
time+=expTime/2.0

#Convert unix time to MJD
JDtime = time/86400.+2440587.5
MJD = JDtime-2400000.5

#convert time to BJD using TEMPO2
parFile = '/home/srmeeker/ARCONS-pipeline/examples/Pal2014_1SWASP_J0002/J0002.par'
BJDdict = processTimestamps(MJD, parFile, workingDir=os.getcwd(), timingProgram='tempo2', bCleanTempDir=True, verbose=False, timingProgramLog='/dev/null')
BJD = BJDdict['baryMjdTimestamps']+2400000.5

#Grab flux arrays from photometry dict
flux=photometryDict['flux']
skyFlux = photometryDict['skyFlux']        

#separate target and ref flux
tar_flux = flux[:,0]
ref_flux = flux[:,1]
flags = photometryDict['flag']

#separate sky fluxes for both
tar_skyFlux = skyFlux[:,0]
ref_skyFlux = skyFlux[:,1]

#subtract sky flux
tar_flux-=tar_skyFlux
ref_flux-=ref_skyFlux

#convert fluxes to counts/sec
tar_flux/=expTime
ref_flux/=expTime

#estimate error as percent std dev of reference star
ref_mean = np.mean(ref_flux)
ref_std = np.std(ref_flux)
percentErr = ref_std/ref_mean

#scale errors to target intensity
tar_errs = percentErr*np.mean(tar_flux)/np.sqrt(tar_flux/ref_mean)
ref_errs = percentErr*ref_flux

#find max and min fluxes with errors for converting to magnitude error bars later
tar_highFlux = tar_flux+tar_errs
tar_lowFlux = tar_flux-tar_errs
ref_highFlux = ref_flux+ref_errs
ref_lowFlux = ref_flux-ref_errs

#plot fluxes in counts/s
plt.errorbar(BJD, tar_flux, yerr = tar_errs, fmt='b.', label='target')
plt.errorbar(BJD, ref_flux, yerr = ref_errs, fmt='g.', label = 'reference')
plt.title(filt)
plt.show()

#convert fluxes to magnitudes
tar_mags = utils.countsToApparentMag(cps=tar_flux,filterName=filt,telescope='Palomar')
ref_mags = utils.countsToApparentMag(cps=ref_flux,filterName=filt,telescope='Palomar')

#get magnitudes of error bars
tar_lowMag = utils.countsToApparentMag(cps=tar_highFlux, filterName=filt,telescope='Palomar')
tar_highMag = utils.countsToApparentMag(cps=tar_lowFlux, filterName=filt,telescope='Palomar')
tar_lowErrs = tar_mags-tar_lowMag
tar_highErrs = tar_highMag-tar_mags

ref_lowMag = utils.countsToApparentMag(cps=ref_highFlux, filterName=filt,telescope='Palomar')
ref_highMag = utils.countsToApparentMag(cps=ref_lowFlux, filterName=filt,telescope='Palomar')
ref_lowErrs = ref_mags-ref_lowMag
ref_highErrs = ref_highMag-ref_mags

#plot light curves in magnitudes
plt.errorbar(BJD,tar_mags,yerr=[tar_lowErrs, tar_highErrs], fmt='b.',label='target')
plt.errorbar(BJD,ref_mags,yerr=[ref_lowErrs, ref_highErrs], fmt='g.',label='ref')
plt.ylim(22,16)
plt.title(filt)
plt.show()

#take out NANs because PHOEBE dies on them
BJD = BJD[np.isnan(tar_mags)==False]
tar_lowErrs = tar_lowErrs[np.isnan(tar_mags)==False]
tar_mags = tar_mags[np.isnan(tar_mags)==False]

BJD = BJD[np.isinf(tar_mags)==False]
tar_lowErrs = tar_lowErrs[np.isinf(tar_mags)==False]
tar_mags = tar_mags[np.isinf(tar_mags)==False]

# output data to PHOEBE formatted txt file
fname = 'PHOEBE_%s_%s.txt'%(objName, filt)
writeToPHOEBEFile(fname, objName, filt, BJD, tar_mags, tar_lowErrs)


