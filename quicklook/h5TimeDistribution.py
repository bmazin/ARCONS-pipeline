
import sys
import time
import struct
import os
from os.path import isfile
#import installed libraries
from numpy import *
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
from tables import *
from matplotlib.backends.backend_pdf import PdfPages

numbins = 2000

datapath = "/Users/ourhero/Documents/Mazinlab/ARCONS-data/LICK2012/20120903/"
file =  "obs_20120903-164458"
datafile = datapath+file

h5file = openFile(str(datafile)+".h5", mode = "r")
print "Loaded file %s.h5"%(datafile)
htable = h5file.root.header.header.read()
print "Exposure time = %i secs"%(htable['exptime'])
bmap = h5file.root.beammap.beamimage.read()
nxpix = shape(bmap)[1]
nypix = shape(bmap)[0]

pp = PdfPages('./%s_timeplots.pdf'%(file))
matplotlib.rcParams['font.size']=6

n=-1

#nxpix = 5
#nypix = 5

allhist = zeros(numbins)

for i in xrange(nypix):
	for j in xrange(nxpix):
			n+=1
			print n
			if bmap[i][j] == '':
				print "skipped pixel at bmap[%i][%i]"%(i,j)
				pass
			else:
				#try:
					photons = h5file.root._f_getChild(bmap[i][j]).read()
					timestamps = array([0])
					#print timestamps
					
					#allphotons= concatenate(h5file.root._f_getChild(bmap[i][j])[:])
					#timestamps = (photons)%pow(2,20)
					for t in xrange(len(photons)):
						#timestamps[t] = (photons)%pow(2,20)
						newtimestamps = array((photons[t]))%pow(2,20)
						newtimestamps += t*1000000
						timestamps = concatenate([timestamps,newtimestamps]) 
					#alltimestamps = (allphotons)%pow(2,20)
					print len(timestamps)
					#subdiffs = zeros(len(alltimestamps-1))
					totaldiffs = zeros(len(timestamps)-1)
					
					#for k in xrange(len(subdiffs)-1):
					#	subdiffs[k] = subdiffs[k] + alltimestamps[k+1] - alltimestamps[k]
					for l in xrange(len(totaldiffs)):
						totaldiffs[l] += timestamps[l+1] - timestamps[l]
					
					print "making histogram for pixel %s"%(bmap[i][j])
					totalhist,totalbins = histogram(totaldiffs,bins=numbins,range=(0,numbins*2.0))
					#subhist,subbins = histogram(subdiffs,bins=numbins,range=(subdiffs.min(),subdiffs.max()))
				
					plt.figure()
					ax = plt.subplot(111)
					ax.set_title(str(bmap[i][j]))

					#plt.bar(totalbins[:-1],totalhist, width = (totalbins[1]-totalbins[0]),bottom=0,label ="With Seconds added")
					plt.plot(totalbins[:-1],totalhist)
					#plt.bar(subbins[:-1],subhist, width = (subbins[1]-subbins[0]),bottom=0,label ="Without Seconds added")
					
					#plt.show()
					#plt.legend()
					pp.savefig()
					
					print "done with pixel %s"%(bmap[i][j])
					allhist += totalhist
				#except:
					#print "failed to make histogram for pixel %s"%(bmap[i][j])
pp.close()


pp = PdfPages('./%s_alltimes.pdf'%(file))
matplotlib.rcParams['font.size']=6
plt.figure()
ax = plt.subplot(111)
ax.set_title("All pixels")
#plt.bar(totalbins[:-1],totalhist, width = (totalbins[1]-totalbins[0]),bottom=0,label ="With Seconds added")
plt.plot(totalbins[:-1],allhist)
print "Finished plot of all pixels stacked"
#plt.legend()
pp.savefig()

pp.close()
h5file.close()
print "DONE"
		