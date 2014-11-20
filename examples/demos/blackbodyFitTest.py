from util import utils
import sys,os
import tables
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from util.ObsFile import ObsFile
from util import MKIDStd
from util.readDict import readDict
from util.rebin import rebin
import matplotlib
from scipy import interpolate
from scipy.optimize.minpack import curve_fit
from numpy import exp

c=3.00E10 #cm/s
h=6.626E-27 #erg*s
k=1.3806488E-16 #erg/K

objectName = 'G158-100'

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(3000,12000)

std = MKIDStd.MKIDStd()
a = std.load(objectName)
a = std.countsToErgs(a)
x = a[:,0]
y = np.array(a[:,1]) #std object spectrum in counts/s/Angs/cm^2
#End MKIDStd loading

plt.plot(x,y*1E15,linewidth=1,color='grey',alpha=0.75)

fraction = 0.65 # set to 1 to fit BB to whole spectrum. 0.2 fits BB to only last 20% of spectrum
fitx = x[(1.0-fraction)*len(x)::]
fity = y[(1.0-fraction)*len(x)::]
plt.plot(fitx, fity*1E15,color='black')

Temp = utils.fitBlackbody(x,y)

newwl = np.arange(4000,20000,100)
Temp, newflux = utils.fitBlackbody(x,y,fraction=fraction,newWvls=newwl)

plt.plot(newwl,newflux*1E15,color = 'blue')
plt.legend(['G158-100 Spectrum','Section used to fit', 'BB Fit'],'upper right', numpoints=1)
plt.xlabel(ur"Wavelength [\AA]")
plt.ylabel(ur"Flux [10$^{-15}$ ergs/s/cm$^{2}$/\AA]")
plt.show()

