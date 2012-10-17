import numpy as np
import scipy as sp
from tables import *
import ConfigParser
import ephem

def arcsec_to_radians(total_arcsec):
    total_degrees = total_arcsec/3600.0
    total_radians = total_degrees*d2r
    return total_radians

def radians_to_arcsec(total_radians):
    total_degrees = total_radians*r2d
    total_arcsec = total_degrees*3600.0
    return total_arcsec
    
d2r = np.pi/180.0
r2d = 180.0/np.pi

# Open the h5 observation file
path = '/Users/kids/desktop/RaDecMap/'
obs_file = 'test_obs.h5'
h5file = openFile(path + obs_file, mode = 'r')

# Set initial variables
ini_file = path + 'LICK2012Initialization.ini'
Config = ConfigParser.ConfigParser()
Config.read(ini_file)
grid_width = Config.getint('ARRAY','GRID_WIDTH')
grid_height = Config.getint('ARRAY','GRID_HEIGHT')
plate_scale = Config.getfloat('ARRAY','PLATE_SCALE')
crpix1 = Config.getfloat('ARRAY','CRPIX1')
crpix2 = Config.getfloat('ARRAY','CRPIX2')
HA_offset = Config.getfloat('OFFSETS','HA_OFFSET')
centroid_RA = Config.get('OFFSETS','CENTROID_RA')
centroid_DEC = Config.get('OFFSETS','CENTROID_DEC')
print 'Grid Width =',grid_width
print 'Grid Height =',grid_height
print 'Plate Scale =',plate_scale
print 'Center of Rotation:',crpix1,',',crpix2
print 'Hour Angle Offset =',HA_offset
print 'Centroid RA =',centroid_RA
print 'Centroid Dec =',centroid_DEC


# Open the centroid positions file
centroids_file = 'centroid_positions.txt'
centroid_x, centroid_y = np.loadtxt(path+centroids_file,unpack='true')

# Create the h5 output file
out_file = 'coords_out.h5'
h5out = openFile(path + out_file, mode = 'w')
ragroup = h5out.createGroup('/','ra', 'RA Map of Array')
decgroup = h5out.createGroup('/','dec', 'DEC Map of Array')
filt1 = Filters(complevel=0, complib='blosc', fletcher32=False)   

# Extract relevant header information from the h5 file
original_lst = h5file.root.header.header.col('lst')[0]
exptime = h5file.root.header.header.col('exptime')[0]
ts = h5file.root.header.header.col('unixtime')[0]
print 'Original LST from telescope:', original_lst
print 'Exptime:', exptime

# Initial RA and LST
centroid_RA_radians = ephem.hours(centroid_RA).real
centroid_RA_arcsec = radians_to_arcsec(centroid_RA_radians)

centroid_DEC_radians = ephem.degrees(centroid_DEC).real
centroid_DEC_arcsec = radians_to_arcsec(centroid_DEC_radians)

original_lst_radians = ephem.hours(original_lst).real
original_lst_seconds = radians_to_arcsec(original_lst_radians)/15.0

rotation_matrix =np.zeros((2,2))
offsets_hypercube = np.zeros(((((grid_height,grid_width,2,exptime)))), dtype = '|S10')
# Calculations done for each second interval
for elapsed_time in range(exptime):
    # Load the image to find the star centroid        
    ra_array = h5out.createCArray(ragroup, 't%i' %elapsed_time, StringAtom(itemsize=10), (grid_height,grid_width), filters = filt1)
    dec_array = h5out.createCArray(decgroup, 't%i' %elapsed_time, StringAtom(itemsize=10), (grid_height,grid_width), filters = filt1)
    # Find an hour angle for each frame, assume center does not move
    current_lst_seconds = original_lst_seconds + elapsed_time
    current_lst_radians = arcsec_to_radians(current_lst_seconds*15.0)
    HA_variable = current_lst_radians - centroid_RA_radians
    HA_static = HA_offset*d2r
    HA_current = HA_variable + HA_static
    # Calculate rotation matrix elements
    rotation_matrix[0][0] = rotation_matrix[1][1] = np.cos(HA_current)
    rotation_matrix[0][1] = -np.sin(HA_current)
    rotation_matrix[1][0] = np.sin(HA_current)

    # Calculate the offsets from center
    x_offsets = np.zeros((grid_height,grid_width))
    y_offsets = np.zeros((grid_height,grid_width))
    rotated_x_offsets = np.zeros((grid_height,grid_width))
    rotated_y_offsets = np.zeros((grid_height,grid_width))
    for y in range(grid_height):
        for x in range(grid_width):
            # Unrotated matrix elements, multiplied by plate scale
            x_offsets[y][x] = plate_scale*(crpix1-x)
            y_offsets[y][x] = plate_scale*(crpix2-y)
            # Apply rotation by hour angle
            rotated_x_offsets[y][x] = centroid_RA_arcsec - plate_scale*(centroid_x[elapsed_time]-crpix1) + rotation_matrix[0][0]*x_offsets[y][x] + rotation_matrix[0][1]*y_offsets[y][x]
            rotated_y_offsets[y][x] = centroid_DEC_arcsec - plate_scale*(centroid_y[elapsed_time]-crpix2) + rotation_matrix[1][0]*x_offsets[y][x] + rotation_matrix[1][1]*y_offsets[y][x]   
            ra_array[y,x] = str(ephem.hours(arcsec_to_radians(rotated_x_offsets[y][x])))
            dec_array[y,x] = str(ephem.degrees(arcsec_to_radians(rotated_y_offsets[y][x])))
            h5out.flush()

