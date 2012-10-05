import numpy as np
import scipy as sp
from tables import *

def DDMMSS_to_seconds(hhmmss_string):
    separator = hhmmss_string.find(':')
    hours = float(hhmmss_string[:separator])
    minutes_seconds = hhmmss_string[(separator+1):]
    separator = minutes_seconds.find(':')
    minutes = float(minutes_seconds[:separator])
    seconds = float(minutes_seconds[(separator+1):])
    total_arcseconds = hours*3600.0 + minutes*60.0 + seconds
    return total_arcseconds

def seconds_to_DDMMSS(arcseconds_float):
    seconds_float = arcseconds_float
    hours = "%02d" %(seconds_float/3600)
    seconds_remainder = seconds_float%3600
    minutes = "%02d" %(seconds_remainder/60)
    seconds = "%04.1f" %(seconds_remainder%60) 
    hhmmss = str(hours) + ':' + str(minutes) + ':' + str(seconds)
    return hhmmss

def seconds_to_radians(total_seconds):
    total_degrees = total_seconds/3600.0
    total_radians = total_degrees*d2r
    return total_radians

def calculate_hourangle(LST_radians,RA_radians):
    HA_radians = LST_radians-RA_radians
    return HA_radians

# Set initial variables
grid_width = 44
grid_height = 46
HA_offset = 19.0
plate_scale = 1.0
centroid_RA = '03:45:26.5'
centroid_DEC = '17:47:53'

# These 2 values should probably be found by finding the center of the gaussian in an image
centroid_x = 22.0
centroid_y = 23.0

d2r = np.pi/180.0
r2d = 180.0/np.pi

# Open the h5 observation file
path = '/Users/kids/desktop/RaDecMap/'
obs_file = 'test_obs.h5'
h5file = openFile(path + obs_file, mode = 'r')

# Extract relevant header information from the h5 file
original_lst = h5file.root.header.header.col('lst')[0]
exptime = h5file.root.header.header.col('exptime')[0]
print 'Original LST from telescope:', original_lst
print 'Exptime:', exptime

# Initial RA and LST
centroid_RA_seconds = DDMMSS_to_seconds(centroid_RA)*15
centroid_RA_radians = seconds_to_radians(centroid_RA_seconds)
centroid_DEC_seconds = DDMMSS_to_seconds(centroid_DEC)
original_lst_seconds = DDMMSS_to_seconds(original_lst)
rotation_matrix =np.zeros((2,2))
offsets_hypercube = np.zeros(((((grid_height,grid_width,2,exptime)))), dtype = '|S10')
# Calculations done for each second interval
for elapsed_time in range(exptime):
    # Find an hour angle for each frame, assume center does not move
    current_lst_seconds = original_lst_seconds + elapsed_time
    current_lst_radians = seconds_to_radians(current_lst_seconds)
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
    hhmmss_RA_offsets = np.zeros(((grid_height,grid_width)),dtype = '|S10')
    for y in range(grid_height):
        for x in range(grid_width):
            # Unrotated matrix elements, multiplied by plate scale
            x_offsets[y][x] = plate_scale*(centroid_x-x)
            y_offsets[y][x] = plate_scale*(centroid_y-y)
            # Apply rotation by hour angle
            rotated_x_offsets[y][x] = rotation_matrix[0][0]*x_offsets[y][x] + rotation_matrix[0][1]*y_offsets[y][x]
            rotated_y_offsets[y][x] = rotation_matrix[1][0]*x_offsets[y][x] + rotation_matrix[1][1]*y_offsets[y][x]
            rotated_x_offsets[y][x] += centroid_RA_seconds
            rotated_y_offsets[y][x] += centroid_DEC_seconds
            offsets_hypercube[y][x][0][elapsed_time] = seconds_to_DDMMSS(rotated_x_offsets[y][x]/15.0)
            offsets_hypercube[y][x][1][elapsed_time] = seconds_to_DDMMSS(rotated_y_offsets[y][x])

