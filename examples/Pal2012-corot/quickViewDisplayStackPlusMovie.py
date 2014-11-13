import numpy as np
import os
import glob
import matplotlib.pylab as plt
from util import utils


files = []
directory = '/Scratch/DisplayStack/PAL2012/Corot18b/ImageStacks/'
fname = "ImageStack*.npz"
outFname = "AllImageStack_20121208_30s_4000-5000_hp.npz"

files = glob.glob(os.path.join(directory, "ImageStack*.npz"))
frames = []
jds = []

for fname in files:
    sDict = np.load(fname)
    #print sDict
    stack = sDict['stack']
    jd = sDict['jd']
    print fname
    print jd
    print np.shape(stack)
    print np.shape(jd)
    #print stack
    #print jd
    for i in xrange(len(stack[0,0,:])):
        print i
        frames.append(stack[:,:,i])
        jds.append(jd[i])
        #plt.matshow(stack[:,:,i])
        #plt.show()

print frames
print jds

np.savez(str(directory+outFname), stack=frames, jd=jds)
utils.makeMovie(frames,frameTitles=jds,cbar=True,outName=directory+'AllImageStack_20121208_30s_4000-5000_hp.gif', normMin=100, normMax=5000)



