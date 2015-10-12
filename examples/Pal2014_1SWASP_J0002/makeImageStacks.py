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

path = '/Scratch/DisplayStack/PAL2014/1SWASP_J0002'
verbose=True
showPlot=True

startWl = 4000
endWl = 5000
hp=True
flat = True
integrationTime=10

outputFN = '%s_%s-%s'%(integrationTime,startWl, endWl)
if flat==True:
    outputFN+='_flat'
if hp==True:
    outputFN+='_hp'
identifier=outputFN

print outputFN

LC=LightCurve(fileID=identifier,path=path,targetName='1SWASP_J0002',run='PAL2014',verbose=True,showPlot=False)
LC.makeImageStack(imageStackFilename='',dt=integrationTime,wvlStart=startWl,wvlStop=endWl,weighted=flat, fluxWeighted=False, getRawCount=False,scaleByEffInt=hp, deadTime=100.e-6)
LC.makeAllCentroidFiles(centroidFilenames=['',''])

'''
#replace everything below here with lightCurve.py example under __main()__
LC = LightCurve(path,verbose=verbose,showPlot=showPlot)
obsFNs = LC.getFNs()
#print 'obs[10]: ',obsFNs[10:11][0].obs()
obsFiles = generateObsObjectList(obsFNs, wvlLowerLimit=startWl, wvlUpperLimit=endWl, loadHotPix=hp, loadWvlCal=True, loadFlatCal=flat, loadSpectralCal=False)
print 'numObs: ',len(obsFiles)
im_dict = getImages(fromObsFile=True, fromPhotonList=False, fromImageStack=False, obsFiles=obsFiles, integrationTime=integrationTime, weighted=flat,  fluxWeighted=False, getRawCount=False)
images=im_dict['images']
pixIntTimes=im_dict['pixIntTimes']
print 'numImages: ',len(images)
    
writeImageStack(images=images, pixIntTimes=pixIntTimes, startTimes=im_dict['startTimes'], intTimes=im_dict['intTimes'], path=path, outputFilename=outputFN)
'''
