import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from util.ObsFile import ObsFile 
from util.FileName import FileName
from util import utils
from util import hgPlot
from cosmic.Cosmic import Cosmic
import tables
from hotpix import hotPixels
import pickle
from interval import interval, inf, imath
import logging, os
import pickle

pickleFile = open("csb2.pkl","r")
fn = pickle.load(pickleFile)
offsets = pickle.load(pickleFile)

fig = plt.figure()
ax = fig.add_subplot(111)

for offset in offsets:
    d = pickle.load(pickleFile)
    print offset, d.keys(),len(d['rows'])
    hg = np.bincount(d['secs'], minlength=300)
    ax.plot(hg,label="offset=%d"%offset)
ax.legend()
plt.title(fn.getComponents())

fig.savefig("csb2plot.png")
pickleFile.close()
