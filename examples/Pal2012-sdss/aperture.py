#!/bin/python

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util.readDict import readDict
from util import utils
import tables
import matplotlib.pyplot as plt
import hotpix.hotPixels as hp
import os
from time import time


def circleAperture(startpx,startpy,radius=3):
    offset = np.floor(radius+0.5)
#    print offset
    xBinsPerPixel = yBinsPerPixel = 50
    totalBinsPerPixel = xBinsPerPixel * yBinsPerPixel
    xPixels = yPixels = 1 + 2*np.floor(0.5 + radius)
#    print xPixels
    xBins = yBins = int(xPixels*xBinsPerPixel)
    pixelsInSquare = int(xPixels * yPixels)
    squareVals = np.linspace(-xPixels/2.,xPixels/2.,xBins)
#    print -xPixels/2,xPixels/2
    squareVals = squareVals[:-1]

    xvals = []
    yvals = []
    for x in range(xBins-1):
        for y in range(yBins-1):
            if squareVals[x]**2 + squareVals[y]**2 <= radius**2:
                xvals.append(squareVals[x])
                yvals.append(squareVals[y])

    counts = np.zeros((int(yPixels),int(xPixels)))
    for x in range(int(xPixels)):
        for y in range(int(yPixels)):
            xPos = x - xPixels/2.0
            yPos = y - yPixels/2.0
            for i in range(len(xvals)):
                if (xvals[i] < xPos + 1.0) and (xvals[i] >= xPos) and (yvals[i] < yPos + 1.0) and (yvals[i] >= yPos):
                    counts[y][x] += 1
    pixelWeights = counts / totalBinsPerPixel
#    print pixelWeights
    pixelMask = np.zeros((46,44))
    for x in range(int(xPixels)):
        for y in range(int(yPixels)):
            if (startpy + y - offset >=0) and (startpx + x - offset >= 0):
                pixelMask[startpy + y - offset][startpx + x - offset] = pixelWeights[y][x]
    return pixelMask
#pixMask=circleAperture(0,0,3)
#print 'pixMask',pixMask
#plt.matshow(pixMask,origin = 'lower')
#plt.show()
