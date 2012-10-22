import numpy as np
from tables import *
import ConfigParser
import ephem

# Config file data, eventually
grid_height = 46
grid_width = 44

# Set the path and h5 file names to be imported
path = '/Users/kids/desktop/CentroidCalc/'
h5beamfile = 'beamimage_sci3gamma.h5'
h5obsfile = 'test_obs.h5'
txtout = 'centroid_positions.txt'

# Make map of where pixels fit in grid in format '/r0/p123/'
h5beam = openFile(path + h5beamfile, mode = 'r')
location_strings = h5beam.root._f_getChild('/beammap/beamimage')

# Open the obs file and extract photon data to calculate com positions in time frames
h5obs = openFile(path + h5obsfile, mode = 'r')
exptime = h5obs.root.header.header.col('exptime')[0]
ts = h5obs.root.header.header.col('unixtime')[0]
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
#            flux_cube[y][x][t] = len(data[t])
            sum_x[t] += x*len(data[t])
            sum_y[t] += y*len(data[t])
            mass[t] += len(data[t])
f = open(path + txtout,'w')
for t in range(exptime):
    com_x[t] = sum_x[t]/mass[t]
    com_y[t] = sum_y[t]/mass[t]
    f.write(str(com_x[t]) + '\t' + str(com_y[t]) + '\n')
f.close()
