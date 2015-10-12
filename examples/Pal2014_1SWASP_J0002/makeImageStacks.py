import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util import utils
from util.readDict import readDict
from scipy import pi
from mpltools	import style
#import figureHeader
import os
import sys
import glob
import tables
from util.getImages import *
from util.popup import *
from photometry.PSFphotometry import PSFphotometry
from photometry.LightCurve import LightCurve

if len(sys.argv)>1:
    filt = sys.argv[1]
else:
    filt = None

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
#LC.makeImageStack(imageStackFilename='',dt=intTime,wvlStart=startWl,wvlStop=endWl,
#                           weighted=True, fluxWeighted=True, getRawCount=False, 
#                           scaleByEffInt=True, deadTime=100.e-6, filterName=filt)
#print "\n\n###########################\nFINISHED MAKING IMAGE STACKS\n#############################\n\n"

#LC.makeAllCentroidFiles(centroidFilenames=['',''])



'''
LC = LightCurve(path,verbose=verbose,showPlot=showPlot)
obsFNs = LC.getFNs()
#print 'obs[10]: ',obsFNs[10:11][0].obs()
obsFiles = generateObsObjectList(obsFNs, wvlLowerLimit=startWl, wvlUpperLimit=endWl, loadHotPix=hp, loadWvlCal=True, loadFlatCal=flat, loadSpectralCal=False)
print 'numObs: ',len(obsFiles)
im_dict = getImages(fromObsFile=True, fromPhotonList=False, fromImageStack=False, obsFiles=obsFiles, integrationTime=10, weighted=flat,  fluxWeighted=False, getRawCount=False)
images=im_dict['images']
pixIntTimes=im_dict['pixIntTimes']
print 'numImages: ',len(images)
    
writeImageStack(images=images, pixIntTimes=pixIntTimes, startTimes=im_dict['startTimes'], intTimes=im_dict['intTimes'], path=path, outputFilename=outputFN)
'''
