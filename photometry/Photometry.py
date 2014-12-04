'''
Author: Alex Walter
Date: Dec 3, 2014

Photometry Super class. Extended by PSFphotometry and AperPhotometry
'''

import numpy as np


class Photometry(object):

    def __init__(self,image,centroid,expTime=None):
        '''

        Inputs:
            image - 2D image of data (0 for dead pixel, shouldn't be any nan's or infs)
                  - Should be fully calibrated, dead time corrected, and scaled up to the effective integration time
            expTime - 2D array of pixel exposure times (0 for dead pixels)
            centroid - list of (col,row) tuples. The first tuple is the target location. The next are reference stars in the field
        '''
        self.image = np.asarray(image)
        self.centroid = centroid
        self.expTime =  np.asarray(expTime) if expTime!=None else 1.0*(np.asarray(image)>0)


















