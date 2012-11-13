'''
This module generates a list of all cal files from the Lick run

Danica Marsden                                           Nov. 1, 2012
'''

#!/bin/env python

import sys, os, shutil

outfile = '/ScienceData/waveCalSolnFiles/calfile_list.txt'
f = open(outfile, 'w')

topdir = '/ScienceData/LICK2012/'
start = 20120901
end = 20120919
n_nights = end - start + 1


for i in range(n_nights):

    dir_end = str(start+i)+'/'

    try:
        night_files = os.listdir(topdir+dir_end)
    except:
        pass
    
    for filename in night_files:
        if 'cal_' in filename:
            f.write(topdir+dir_end+filename+'\n')

f.close()












