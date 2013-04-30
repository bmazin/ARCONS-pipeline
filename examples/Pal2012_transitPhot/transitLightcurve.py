'''Julian's version of sdsslightcurve.py'''

import numpy as np
import matplotlib.pyplot as plt
from util import utils
from util.readDict import readDict

param = readDict()
param.read_from_file('photometryParams.dict')

fileName = param['npzLoadFitFile']
framesPerFile = param['FramesPerFile']
print 'Frames per file: ', framesPerFile
#totalNumFiles = param['TotalNumFiles']
#numFrames = framesPerFile*TotalNumFiles
intTime = param['integrationTime']

t = np.load(fileName)
params = t['params']
jd = t['jd']
t.close()
heights = params[:,1]
widths = params[:,4]
xpos = params[:,2]
ypos = params[:,3]

fig = plt.figure()
ax = fig.add_subplot(111)
curve = heights*widths**2
curve /= np.median(curve)

fwhm = 2*np.sqrt(2*np.log(2))*widths#pixels

#jd = jd[curve<=500]
#curve = curve[curve<=500]
#jd = jd[curve>=100]
#curve = curve[curve>=100]

#curve/=np.median(curve)
#fwhm/=np.median(fwhm)
ax.plot(jd,curve,'k.')
ax.set_title(fileName)
plt.show()
