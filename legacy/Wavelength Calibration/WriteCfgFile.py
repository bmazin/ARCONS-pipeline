#####
#
#WriteCfgFile.py
#Takes list of resonators from CalFitOut and makes new cfg file with whichever you choose (flagged 0 or 1, etc)
#
#####

import sys
import mpfit
import numpy as np
import matplotlib
#matplotlib.use('PDF')
import matplotlib.pyplot as plt
from tables import *
from matplotlib.backends.backend_pdf import PdfPages

path = './223423/'

infile = path + 'CalFitOut.txt'

pixels = []
m = []
b = []
flags=[]

f = open(infile, 'r')
for line in f:
		pixel,min,bin,flag = line.split('\t')
		pn, ts = pixel.split('t')
		min = float(min.strip())
		bin = float(bin.strip())
		flag = int(flag.strip())

		pixels.append(pn)
		m.append(min)
		b.append(bin)
		flags.append(flag)
		
		print pn
		print min
		print bin
		print repr(flag)
f.close()

#fitting will output text file with path and file, and list of pixels and fit parameters
#read in text list to arrays for pixels and fit parameters

'''
#temporary solution: write out pixel arrays manually here
pixels.append('r1/p5/')
pixels.append('r1/p11/')
pixels.append('r1/p15/')
pixels.append('r1/p24/')
pixels.append('r1/p27/')
pixels.append('r1/p38/')
pixels.append('r1/p95/')
pixels.append('r1/p124/')
pixels.append('r1/p125/')
pixels.append('r1/p145/')
pixels.append('r1/p160/')
pixels.append('r1/p162/')
pixels.append('r1/p169/')
pixels.append('r1/p174/')

m.append(-0.0046539)
m.append(-0.0056801)
m.append(-0.004788)
m.append(-0.0053794)
m.append(-0.0043035)
m.append(-0.0041663)
m.append(-0.0045049)
m.append(-0.0048132)
m.append(-0.0049837)
m.append(-0.0049299)
m.append(-0.0045272)
m.append(-0.0046421)
m.append(-0.0049432)
m.append(-0.0062)

b.append(9.0684)
b.append(11.327)
b.append(9.6615)
b.append(10.493)
b.append(8.3566)
b.append(8.0057)
b.append(9.0936)
b.append(8.8634)
b.append(9.8249)
b.append(9.8066)
b.append(8.6977)
b.append(9.601)
b.append(9.8396)
b.append(12.6934)
'''

outfile = path + 'badpix.cfg'

fout = open(outfile, mode = "w")

count = 0

for k in range(len(pixels)):
	if flags[k]==0:
					print "Writing pixel ", pixels[k]
					fout.write("pixel="+pixels[k]+"\n")
					count +=1
print "Wrote  ", count, "pixels."
fout.write("dir=/Users/ourhero/Documents/python/MazinLab/Arcons/\n")
fout.write("20110804/obs_20110804-223726.h5")
fout.close()