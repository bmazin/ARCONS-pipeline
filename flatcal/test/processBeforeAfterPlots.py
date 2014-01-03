import numpy as np
import matplotlib.pyplot as plt
from util.popup import PopUp,plotArray
import itertools


#from http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

#from http://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z


imgDict = np.load('beforeAfterImgs.npz')
beforeImg = imgDict['beforeImg']
afterImg = imgDict['afterImg']
deadBeforeImg = beforeImg == 0
deadAfterImg = afterImg == 0


nRows,nCols = np.shape(beforeImg)

#beforeList = beforeImg[beforeImg != 0]
#
#hotSigmaCutoff = 2.
#
#hotBefore = beforeImg > np.mean(beforeList)+hotSigmaCutoff*np.std(beforeList)
#beforeImg[hotBefore] = 0.
#afterImg[hotBefore] = 0.

beforeList = beforeImg[beforeImg != 0]
afterList = afterImg[afterImg != 0]


plotArray(title='without flatcal',image=beforeImg)
plotArray(title='with flatcal',image=afterImg)

beforeHist,beforeHistEdges = np.histogram(beforeList,bins=200)
afterHist,afterHistEdges = np.histogram(afterList,bins=300)


print 'before:'
print 'count',len(beforeList)
print 'sdev',np.std(beforeList)

print 'after:'
print 'count',len(afterList)
print 'sdev',np.std(afterList)

xx,yy = np.meshgrid(np.arange(nCols),np.arange(nRows))
print xx

z = beforeImg.ravel()
x = xx.ravel()
y = yy.ravel()

x = x[z != 0]
y = y[z != 0]
z = z[z != 0]


# Fit a 3rd order, 2d polynomial
m = polyfit2d(x,y,z)
print m

# Evaluate it on a grid...
beforeImgFit = polyval2d(xx, yy, m)
plotArray(beforeImgFit,vmin=0,title='poly fit to non-flatcal image')

afterImgSub = afterImg - np.mean(afterList)
afterImgSub[deadAfterImg] = 0

beforeImgSub = beforeImg - beforeImgFit
beforeImgSub[deadBeforeImg] = 0

plotArray(beforeImgSub,title='without flatcal minus poly fit')
plotArray(afterImgSub,title='with flatcal minus mean')

afterImgSub[deadAfterImg] = np.nan
beforeImgSub[deadBeforeImg] = np.nan

subAfterList = afterImgSub[~np.isnan(afterImgSub)]
subBeforeList = beforeImgSub[~np.isnan(beforeImgSub)]
print np.shape(subBeforeList),np.shape(subAfterList)
subBeforeHist,subBeforeHistEdges = np.histogram(subBeforeList,bins=300)
subAfterHist,subAfterHistEdges = np.histogram(subAfterList,bins=300)

# Plot
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(beforeHistEdges[0:-1],beforeHist,label='without flatcal')
ax.plot(afterHistEdges[0:-1],afterHist,label='with flatcal')
ax.set_title('Distribution of pixel counts')
ax.set_xlabel('Counts')
ax.set_ylabel('Num of Pixels')
ax.legend()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title('Distribution of pixel counts after subtracting')
ax.plot(subBeforeHistEdges[0:-1],subBeforeHist,label='non-flatcal image - poly fit')
ax.plot(subAfterHistEdges[0:-1],subAfterHist,label='flatcal - mean count')
ax.set_xlabel('Counts')
ax.set_ylabel('Num of Pixels')
ax.legend()
plt.show()

