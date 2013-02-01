import numpy as np
from tables import *
import sys
import ephem
import PyGuide as pg
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import os
from PyQt4.QtGui import *

class MouseMonitor():
    def __init__(self):
        pass

    def on_click(self,event):
	self.xyguess = [event.xdata,event.ydata]
   
    def connect(self):
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
	

# Set the path and h5 file names to be imported
h5beamfile = '/ScienceData/sci4alpha/sci4_beammap_palomar.h5'

app = QApplication(sys.argv)
startingSaveDir = os.getcwd()
savePath = QFileDialog.getExistingDirectory(None,'Choose Save Path',startingSaveDir, QFileDialog.ShowDirsOnly)
startingh5Dir = '/ScienceData/PAL2012/'
h5obsfile =str(QFileDialog.getOpenFileName(None, 'Choose Observation File',startingh5Dir, filter=str("H5 (*.h5)")))
n0 = h5obsfile.rfind('/')
n1 = h5obsfile.find('.')
centroid_list_identifier = h5obsfile[n0+1:n1]
txtout = 'centroids_' + centroid_list_identifier + '.txt'

grid_width = 44
grid_height = 46

# Make map of where pixels fit in grid in format '/r0/p123/'
h5beam = openFile(h5beamfile, mode = 'r')
location_strings = h5beam.root._f_getChild('/beammap/beamimage')

# Open the obs file and extract photon data to calculate com positions in time frames
h5obs = openFile(h5obsfile, mode = 'r')
exptime = h5obs.root.header.header.col('exptime')[0]

try:
    ts = h5obs.root.header.header.col('unixtime')[0]
except KeyError:
    print 'Using "ut" instead of "unixtime" in header'
    ts = h5obs.root.header.header.col('ut')[0]

flux_cube = np.zeros(((grid_height,grid_width,exptime)))
for y in range(grid_height):
    for x in range(grid_width):
        pn = location_strings[y][x] + 't' + str(int(ts))
        data = h5obs.root._f_getChild(pn).read()
        for t in range(exptime):
            flux_cube[y][x][t] = len(data[t])

# Specify pixel info CCDInfo(bias,readNoise,ccdGain,satLevel)
ccd = pg.CCDInfo(0,0.00001,1,2500)

# Create a mask of bad pixels, 0 for valid, 1 for invalid
mask = np.zeros((grid_height,grid_width))

# Create a mask for saturated pixels, 0 for okay, 1 for saturated
satMask = np.zeros((grid_height,grid_width))

integration_time= 10
xyguess=[0,0]

f = open(savePath+ '/' + txtout,'w')
for t in range(int(exptime/integration_time)):
    map = MouseMonitor()
    map.fig = plt.figure()
    pltmat = np.zeros((grid_height,grid_width))
    for y in range(grid_height):
        for x in range(grid_width):
	    for i in range(integration_time):
                pltmat[y][x]+=flux_cube[y][x][int(t*integration_time+i)]/integration_time
    map.ax = map.fig.add_subplot(111)
    map.ax.matshow(pltmat,cmap = plt.cm.gray, origin = 'lower')
    map.connect()
    plt.show()
    try:
    	xyguess = map.xyguess
    except AttributeError:
	pass
    print 'Guess = ' + str(xyguess)
    pyguide_output = pg.centroid(pltmat,mask,mask,xyguess,10,ccd,0,False)
    try:
        xycenter = pyguide_output.xyCtr
        print 'Calculated = ' + str(xycenter)
        f=open(savePath+ '/' + txtout,'a')
        f.write(str(xycenter[0]) + '\t' + str(xycenter[1]) + '\n')
        f.close()
    except TypeError:
        print 'Cannot centroid, using guess'
        f=open(savePath+ '/' + txtout,'a')
        f.write(str(xyguess[0]) + '\t' + str(xyguess[1]) + '\n')
        f.close()

