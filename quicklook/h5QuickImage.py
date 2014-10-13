#!/bin/python

import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from util.ObsFile import ObsFile
from PIL import Image
from PIL import ImageDraw

color = True

startwl = 3000
stopwl = 10000
starttime=0
endtime=-1

#datapath = "/Users/ourhero/Documents/Mazinlab/ARCONS-data/PAL2012/20121206/"
#datafiles = ["obs_20121207-052426.h5","obs_20121207-052928.h5"]

datapath = "/Users/ourhero/Documents/Mazinlab/ARCONS-data/PAL2012/20121207/"
datafiles = ["obs_20121208-025525.h5"]

calpath =  "/Users/ourhero/Documents/Mazinlab/ARCONS-data/PAL2012/waveCalSolnFiles/"
calsolnfile = "calsol_20121207-023127.h5"
calfile = calpath + calsolnfile

n=0

for datafile in datafiles:
    n+=1
    print "imaging file ",n, " out of ",len(datafiles)
    obsfile = datapath + datafile
    print "reading in file ", obsfile
    print "using cal file ", calfile

    h5 = ObsFile(obsfile)
    print "ObsFile loaded"
    h5.loadWvlCalFile(calfile)
    print "WvlCal loaded"

    nypix = h5.nRow
    nxpix = h5.nCol

    image = np.empty((nypix,nxpix),dtype=int)
    colors = np.empty((nypix,nxpix),dtype = int)

    print "Making image"

    for i in xrange(nypix):
        for j in xrange(nxpix):
            photons = h5.getPixelWvlList(i,j,starttime,endtime)
            ind = np.where(np.logical_and((photons > startwl),(photons < stopwl)))
            #print "uncut counts = ", len(photons)
            photons = photons[ind]
            image[i,j] = len(photons)
            #print "cut counts = ", image[i,j]
            if color == True:
                try:
                    colors[i,j] = np.median(photons)
                except:
                    colors[i,j] = 10000
                    
        print "..." 
    print "Finished constructing image"
    if n==1:
        totalimage = image
        if color == True:
            totalcolors = colors
    else:
        totalimage += image
        if color == True:
            totalcolors = (n*totalcolors + colors)/(n+1)
    h5.__del__()


print "subtracting off median counts of ", np.median(totalimage[totalimage>0])

totalimage - np.median(totalimage[totalimage>0])
#vmin = np.mean(totalimage)-1.0*np.std(totalimage)
vmin = 0
vmax = np.mean(totalimage)+2.0*np.std(totalimage)

im = Image.new("RGB", (nxpix,nypix))
draw = ImageDraw.ImageDraw(im)

if color == True:
    totalimage[totalimage<vmin] = vmin
    totalimage[totalimage>vmax] = vmax
    totalimage -= vmin
    V = (totalimage/vmax)*80
    #S = np.ones_like(V)
    #S*=0.25 # turn down saturation so colors do not dominate
    H = np.ones_like(V)
    for y in xrange(nypix):
        for x in xrange(nxpix):
            if colors[y][x] < 3900:
                H[y][x] = 300
            elif colors[y][x] > 7500:
                H[y][x] = 0
            else:
                H[y][x] = colors[y][x] * -0.0732 + 556.1
            h = H[y][x]
            s = V[y][x] #scale S with V so brightest pixels show most color
            v = V[y][x]
            colorstring = "hsl(%i,%i%%,%i%%)"%(h,s,v)
            draw.point((x,nypix-y),colorstring)

    im.show()

#    HSV = np.dstack((H,S,V))
#    RGB = hsv_to_rgb(HSV)   
#    plt.imshow(RGB,interpolation = 'nearest')

else:
    RGB = totalimage
    plt.gray()
    plt.matshow(RGB,vmin = vmin, vmax = vmax)
    plt.show()
