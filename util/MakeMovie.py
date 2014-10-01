from ObsFile import ObsFile
import numpy as np
import tables
import sys,os
import matplotlib.pyplot as plt
from images2gif import writeGif
import time

'''
Given an h5 observation file, makes images of every 1s frame and saves them to .png format,
to be used later to make into a movie (with iMovie).  GIF functionality not implemented compeltely.
'''

gif = False

manual = True

vmin = 0
vmax = 50

#change to take filename as argument from command line
datapath = "/Users/ourhero/Documents/Mazinlab/ARCONS-data/PAL2012/"

frame = 0.1 #s, only supports values less than 1s

file = "obs_20121203-070614.h5"
base = file.split(".")[0]

filename = datapath + file
outpath = datapath + "images/"+ base + "/"

if os.path.exists(outpath) == False:
    os.mkdir(outpath)
    print "outpath not found, creating ", outpath

print "opening file ", filename

h5f = tables.openFile(filename,mode='r')
bmap = h5f.root.beammap.beamimage.read()
nypix = bmap.shape[0]
nxpix = bmap.shape[1]
header = h5f.root.header.header.read()
t = int(header["exptime"])
nframes=t/frame
fps = 1.0/frame

print "obs time = ", t
print "number of frames = ", nframes

img = np.zeros((nframes,nypix,nxpix))
fig = plt.figure()
plt.gray()

for i in xrange(t):
    print "displaying second ", i
    for n in xrange(int(fps)):
        fnum = n+fps*i
        print "displaying frame ", fnum
        timestart = time.time()
        for y in xrange(nypix):
            for x in xrange(nxpix):                
                packets = h5f.root._f_getChild(bmap[y][x])[i]
                times = np.bitwise_and(packets, int(20*'1',2))
                ind = np.where(times > n*frame*1E6)
                times = times[ind]
                ind = np.where(times < (n+1)*frame*1E6)
                times = times[ind]
                img[fnum,y,x] = len(times)
        ax = fig.add_subplot(111)
        if manual == True:
            ax.matshow(img[fnum],vmin = vmin, vmax = vmax)
        else:
            ax.matshow(img[fnum],vmax = np.mean(img[fnum])+2*np.std(img[fnum]))
        fig.savefig(outpath + base+"_%i.png"%(fnum))

        del ax
        fig.clf()

        timefinish = time.time()
        print "runtime = ", timefinish-timestart

plt.close(fig)
h5f.close()
if gif == True:
    print "Turning file into GIF..."
    writeGif(outpath + base + ".gif",img)

print "DONE"
