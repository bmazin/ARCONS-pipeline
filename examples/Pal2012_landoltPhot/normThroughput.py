import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

#normalizes a throughput file with the correct curve shape to the correct absolute height as measured by 
#a file that was generated on a better night of data.

maxfile = '/home/srmeeker/scratch/standards/sdss j0926_throughput.npz'

curvefile = '/home/srmeeker/scratch/standards/G158-100_throughput.npz'

maxDict = np.load(maxfile)
curveDict = np.load(curvefile)

maxThru = maxDict['throughput']
curveThru = curveDict['throughput']
wvls = maxDict['wvls']

print wvls
print maxThru
print curveThru


peak = max(maxThru[(wvls>3900) & (wvls<8000)])
peakwvl = wvls[maxThru==max(maxThru[(wvls>3900) & (wvls<8000)])]

curveThru /= curveThru[wvls == peakwvl]
print curveThru

curveThru *= peak
print curveThru

plt.plot(wvls,curveThru)
plt.xlim(3900,13000)
plt.ylim(0,0.1)
plt.show()

outarr = np.empty((len(wvls),2),dtype=float)
outarr[:,0]=wvls
outarr[:,1]=curveThru
#save throughput curve to file
np.savetxt("throughput.txt", outarr)
