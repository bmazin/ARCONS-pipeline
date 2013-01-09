import numpy as np
from tables import *
import ConfigParser
import ephem
import PyGuide as pg
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Initial guess of centroid location
xguess = 14.0
yguess = 14.0
xyguess = [yguess,xguess]

# Set the path and h5 file names to be imported
path = '/home/pszypryt/Scratchwork/'
h5beamfile = 'beamimage_sci3gamma.h5'
h5obsfile = 'obs_20120904-230531.h5'
n = h5obsfile.find('.')
centroid_list_identifier = h5obsfile[0:n]
txtout = 'centroids_' + centroid_list_identifier + '.txt'
ini_file = path + 'Lick2012Initialization.ini'

# Config File Information
Config = ConfigParser.ConfigParser()
Config.read(ini_file)
grid_width = Config.getint('ARRAY','GRID_WIDTH')

grid_height = Config.getint('ARRAY','GRID_HEIGHT')

# Make map of where pixels fit in grid in format '/r0/p123/'
h5beam = openFile(path + h5beamfile, mode = 'r')
location_strings = h5beam.root._f_getChild('/beammap/beamimage')

# Open the obs file and extract photon data to calculate com positions in time frames
h5obs = openFile(path + h5obsfile, mode = 'r')
exptime = h5obs.root.header.header.col('exptime')[0]
try:
    ts = h5obs.root.header.header.col('unixtime')[0]
except KeyError:
    print 'Using "ut" instead of "unixtime" in header'
    ts = h5obs.root.header.header.col('ut')[0]    
flux_cube = np.zeros(((grid_height,grid_width,exptime)))
sum_x = np.zeros(exptime)
sum_y = np.zeros(exptime)
mass = np.zeros(exptime)
com_x = np.zeros(exptime)
com_y = np.zeros(exptime)
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

f = open(path + txtout,'w')
for t in range(exptime):
    grid_to_centroid = np.zeros((grid_height,grid_width))
    for y in range(grid_height):
        for x in range(grid_width):
            grid_to_centroid[y][x] = flux_cube[y][x][t]          
    pyguide_output = pg.centroid(grid_to_centroid,mask,mask,xyguess,10,ccd,0,False)
    xycenter = pyguide_output.xyCtr
    f.write(str(xycenter[0]) + '\t' + str(xycenter[1]) + '\n')
    if (t==20):
        plt.imshow(grid_to_centroid,cmap = plt.cm.gray)
f.close()
plt.show()
