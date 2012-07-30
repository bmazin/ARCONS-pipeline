import sys
import mpfit
import numpy as np
import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
from tables import *
from matplotlib.backends.backend_pdf import PdfPages

path = './045533/'
datapath = './20110804/'

infile = datapath + 'obs_20110805-045533.h5'
outfile = datapath + 'cal_20110805-045533.h5'

pixels = []
m = []
b = []
flags=[]

mfull = []
bfull = []
flagsfull=[]

fbasic = open(path+"CalFitOut_all_basic.txt", 'r')
ffull = open(path+"CalFitOut_all_full.txt", 'r')

for line in fbasic:
		pixel,min,bin,flag = line.split('\t')
		pn, ts = pixel.split('t')
		min = float(min.strip())
		bin = float(bin.strip())
		flag = int(flag.strip())

		pixels.append(pn)
		m.append(min)
		b.append(bin)
		flags.append(flag)
		
for line in ffull:
		pixel,min,bin,flag = line.split('\t')
		pn, ts = pixel.split('t')
		minfull = float(min.strip())
		binfull = float(bin.strip())
		flagfull = int(flag.strip())
		
		#pixels.append(pn)
		mfull.append(minfull)
		bfull.append(binfull)
		flagsfull.append(flagfull)

ffull.close()
fbasic.close()

#fitting will output text file with path and file, and list of pixels and fit parameters
#read in text list to arrays for pixels and fit parameters

h5in = openFile(infile,mode = "r")
bmap = h5in.root.beammap.beamimage.read()



h5out = openFile(outfile, mode = "w")
bgroup = h5out.createGroup('/','calparams','Table of calibration parameters for each pixel')

filt1 = Filters(complevel=1, complib='blosc', fletcher32=False)   # without minimal compression the files sizes are ridiculous...

h5out.copyNode(h5in.root.beammap,newparent=h5out.root,recursive=True,overwrite=True)
h5in.close()
print "Copied over beammap"

calarray = h5out.createCArray(bgroup, 'energyscale', Float32Atom(), (32,32,2), filters=filt1)
flagarray = h5out.createCArray(bgroup, 'qualityflag', Int32Atom(), (32,32), filters=filt1)

countfull = 0
countbasic = 0

for k in range(len(pixels)):
	if flags[k]==1:
		for i in range(32):
			for j in range(32):
				pix = bmap[i][j].split('/')
				pixtest = '/'+str(pix[1]+'/'+pix[2]+'/')
				if pixels[k] == pixtest:
					print "Writing pixel ", pixels[k]
					outi = i
					outj = j
					print "i,j = ", i, j
					calarray[outi,outj] = np.array([m[k],b[k]])
					flagarray[outi,outj] = 1
					h5out.flush()
					countbasic +=1
	if flagsfull[k]==1:
		for i in range(32):
			for j in range(32):
				pix = bmap[i][j].split('/')
				pixtest = '/'+str(pix[1]+'/'+pix[2]+'/')
				#print pixels[k]
				#print pixtest
				if pixels[k] == pixtest:
					#print "Writing pixel ", pixels[k]
					outi = i
					outj = j
					#print "i,j = ", i, j
					calarray[outi,outj] = np.array([mfull[k],bfull[k]])
					flagarray[outi,outj] = 2
					h5out.flush()
					countfull +=1
					
print "Wrote basic cal data for ", countbasic, " pixels."
print "Wrote full cal data for ", countfull, " pixels."
h5out.close()