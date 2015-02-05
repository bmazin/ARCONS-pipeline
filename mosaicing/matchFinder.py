import numpy as np
import tables
import sys
import ephem
import PyGuide as pg
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util.ObsFile import ObsFile 
import os
from PyQt4 import QtGui
from PyQt4 import QtCore
from PyQt4.QtGui import *
import hotpix.hotPixels as hp
from tables import *
from util.FileName import FileName
from photometry.PSFphotometry import PSFphotometry
from util import utils
from util.popup import PopUp
import astrometry.CentroidCalc as cc
import util.ObsFileSeq as ofs

import time


#added by Neil vvvvvv 2/3/2015

class MouseMonitor():
    def __init__(self):
        pass

    def on_click(self,event):
        if event.inaxes is self.ax1:
            self.xyguess1 = [event.xdata,event.ydata]
            print 'Clicked: ',self.xyguess1
        
        elif event.inaxes is self.ax2:
            self.xyguess2 = [event.xdata,event.ydata]
            print 'Clicked: ',self.xyguess2
    
    def on_scroll_cbar(self,event):
        if event.inaxes is self.fig1.cbar.ax:
            increment=0.05
            currentClim = self.fig1.cbar.mappable.get_clim()
            currentRange = currentClim[1]-currentClim[0]
            if event.button == 'up':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]+increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]+increment*currentRange)
            if event.button == 'down':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]-increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]-increment*currentRange)
            self.fig1.cbar.mappable.set_clim(newClim)
            self.fig1.canvas.draw()
           
        elif event.inaxes is self.fig2.cbar.ax:
            increment=0.05
            currentClim = self.fig2.cbar.mappable.get_clim()
            currentRange = currentClim[1]-currentClim[0]
            if event.button == 'up':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]+increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]+increment*currentRange)
            if event.button == 'down':
                if QtGui.QApplication.keyboardModifiers()==QtCore.Qt.ControlModifier:
                    newClim = (currentClim[0]-increment*currentRange,currentClim[1])
                elif QtGui.QApplication.keyboardModifiers()==QtCore.Qt.NoModifier:
                    newClim = (currentClim[0],currentClim[1]-increment*currentRange)
            self.fig2.cbar.mappable.set_clim(newClim)
            self.fig2.canvas.draw()
    
    def connect(self):
        self.cid1 = self.fig1.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig1.cbar = self.fig1.colorbar(self.handleMatshow1)
        cid1 = self.fig1.canvas.mpl_connect('scroll_event', self.on_scroll_cbar)

        self.cid2 = self.fig2.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig2.cbar = self.fig2.colorbar(self.handleMatshow2)
        cid2 = self.fig2.canvas.mpl_connect('scroll_event', self.on_scroll_cbar)


def getUserObjectGuess(images, norm = None):
    '''
    This is designed to allow the user to look at two frames, determine whether they contain matching stars,
    and if so click on their position withing the image.
    '''
    flagList = np.array([1,1])
    xyguess1 = [0,0]
    xyguess2 = [0,0]

    image1 = images[0]
    image2 = images[1]

    map = MouseMonitor()
    map.fig1 = plt.figure(1)
    map.ax1 = map.fig1.add_subplot(111)
    map.ax1.set_title('Star Position Guess')
    map.handleMatshow1 = map.ax1.matshow(image1,cmap = mpl.cm.gnuplot2, origin = 'lower', norm=norm)
    
    map.fig2 = plt.figure(2)
    map.ax2 = map.fig2.add_subplot(111)
    map.ax2.set_title('Star Position Guess')
    map.handleMatshow2 = map.ax2.matshow(image2,cmap = mpl.cm.gnuplot2, origin = 'lower', norm=norm)

    map.connect()
    
    plt.show() 
    
    try:
        xyguess1 = map.xyguess1
        print 'Guess1 = ' + str(xyguess1)
        flagList[0]=0
    except AttributeError:
        pass

    try:     
        xyguess2 = map.xyguess2
        print 'Guess2 = ' + str(xyguess2)
        flagList[1]=0
    except AttributeError:
        pass
    
    xyguesses = np.array([xyguess1,xyguess2])    

    return xyguesses, flagList

def ObjectFinder(frames, data, RA = None, Dec = None, radiusOfSearch=10, usePsfFit=False):
    '''
    This allows the user to determine the exact positioning of a star withing a given frame, given that
    the star exists in two frames. two images pop up, if they contain the same star, click on it in both images,
    and the exact pixel positioning will be determine via centroiding. Some parts are fudged for now.

    Inputs:
        frames - This is an array of images corresponding to the various frames used for the mosaic. note that 
                 the way this is written as of now requires that the first  and last frame include the same star. 
                 I think this is typical so it should not be a problem.

        data - this is an array of dictionaries corresponding to each frame. This is generated from
               getFrameDict() which I have included in the ObsFileSeq class. This is required.

        RA/Dec - these are not used as of now.
    
    Outputs:
        
        all the parameters that are needed for setRm()

    Ben also suggested that a cross-correlating technique may be useful for matching frames with no star in them.
    What do you guys think? I need to look into it more but I believe that this code can be expanded to also do
    cross-correlation. - Neil
    '''

    offsRA = []
    offsDec = []
    iframe = [] 
    meanTime = []  
    
    for i in range(len(data)):
        offsRA.append(data[i]["offsRA"])
        offsDec.append(data[i]["offsDec"])
        iframe.append(data[i]["iframe"])
        meanTime.append(data[i]["meanTime"])

    offsRA = np.array(offsRA)
    offsDec = np.array(offsDec)
    iframe = np.array(iframe)
    meanTime = np.array(meanTime)
    
    dpp = []
    ang = []
    
    for i in range(1, len(frames)):
        images = np.array([frames[i-1], frames[i]])    
        oRA = np.array([offsRA[i-1], offsRA[i]])
        oDec = np.array([offsDec[i-1], offsDec[i]])
        print 'Looking for matching Stars...'

        print 'frame: ', iframe[i-1], 'Offset RA: ', offsRA[i-1], 'Offset Dec: ', offsDec[i-1]
        print 'frame: ', iframe[i], 'Offset RA: ', offsRA[i], 'Offset Dec: ', offsDec[i]         
        
        xyguesses, flagList = getUserObjectGuess(images)

        if flagList[1]==0 and flagList[0]==0:
            print 'match found! - determining centroid positions'
                
            xycenter1, flag1 = cc.centroidImage(images[1], xyguesses[1], radiusOfSearch=radiusOfSearch, doDS9=False, usePsfFit=usePsfFit)
            xycenter0, flag0 = cc.centroidImage(images[0], xyguesses[0], radiusOfSearch=radiusOfSearch, doDS9=False, usePsfFit=usePsfFit)
            print 'Success! Matching stars at: ', xycenter0, xycenter1                
                
            rc1 = np.array(xycenter1)
            rc0 = np.array(xycenter0)
            dCol = rc1[0] - rc0[0]
            dRow = rc1[1] - rc0[1]                
            dPix = math.sqrt(((rc1-rc0)**2).sum())

            #center ra,dec of fram calculated from offsets
            dRA = oRA[1]-oRA[0] #arcseconds
            dDec = oDec[1]-oDec[0]
            dDeg = math.sqrt(dRA**2+dDec**2)/3600 #degrees
                
            degPerPix = dDeg/dPix #plate scale
    
            #rotation
            thetaPix = math.atan2(dCol,dRow) #angle from verticle
            thetaSky = math.atan2(dRA, dDec) #angle from north
            theta = thetaPix-thetaSky        #degrees
                
            dpp.append(degPerPix)
            ang.append(theta)
        
        elif flagList[1]==1 or flagList[0]==1:
            print 'no star found' 

    dpp = np.array(dpp)
    #print dpp
    degPerPix = np.mean(dpp)
    #print degPerPix
    ang = np.array(ang)
    #print ang
    theta = np.mean(ang)  
    #print theta  
     
    ## Pick two frames where the ra,dec offset is zero,
    # usually the beginning and ending frames
    
    print 'Matching stars from the first and last frames'    
    
    images = [frames[0], frames[-1]]

    print 'frame: ', iframe[0], 'Offset RA: ', offsRA[0], 'Offset Dec: ', offsDec[0]
    print 'frame: ', iframe[-1], 'Offset RA: ', offsRA[-1], 'Offset Dec: ', offsDec[-1]  
    
    xyguesses, flagList = getUserObjectGuess(images)
    
    xycenter1, flag1 = cc.centroidImage(images[1], xyguesses[1], radiusOfSearch=radiusOfSearch, doDS9=False, usePsfFit=usePsfFit)
    xycenter0, flag0 = cc.centroidImage(images[0], xyguesses[0], radiusOfSearch=radiusOfSearch, doDS9=False, usePsfFit=usePsfFit)
    
    print 'Success! Matching stars at: ', xycenter0, xycenter1        

    rcA = np.array(xycenter1)
    rcB = np.array(xycenter0)
    sct = math.cos(theta)*degPerPix    
    sst = math.sin(theta)*degPerPix
    # This rotation matrix converts from row,col to ra,dec in degrees
    rm = np.array([[sct,-sst],[sst,sct]])
    rdA = rm.dot(rcA)
    rdB = rm.dot(rcB)
    deltaRa = rdB[0]-rdA[0]
    deltaTime = meanTime[-1] - meanTime[0]
    raArcsecPerSec = 3600*deltaRa/deltaTime
    
    return degPerPix, theta, raArcsecPerSec

