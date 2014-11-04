'''
Author: Danica Marsden                          Date: June 4, 2013

'''

import sys, os
import time
import tables
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import mpfit
import smooth
from utils import *
from tables import *
from fitFunctions import *
from matplotlib.backends.backend_pdf import PdfPages


#ctsR = open('countRateResolutionData2.txt', 'r')
ctsR = open('countRateResolutionData.txt', 'r')

arcons_cts = []
r_avg = []
for line in ctsR:
    if '#' in line:
        continue

    wholeline = line.strip()
    print wholeline.split(' ')
    arcons_cts.append(float(wholeline.split(' ')[1]))
    r_avg.append(float(wholeline.split(' ')[3]))
    
arcons_cts = np.array(arcons_cts)
r_avg = np.array(r_avg)
ctsR.close()

polycoeffs = np.polyfit(arcons_cts, r_avg, 2)
print polycoeffs
##[  7.76649502e-07  -3.97849064e-03   8.58104476e+00]
r_fit = np.polyval(polycoeffs, arcons_cts)

fig = plt.figure()
plt.plot(arcons_cts, r_avg)
plt.plot(arcons_cts, r_fit, 'r')
plt.show()
plt.clf()






