import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from sdssgaussfitter import twodgaussian


#rawfile = '/home/srmeeker/scratch/standards/Landolt9542_raw_0.npz'

#rawfile = '/home/srmeeker/scratch/standards/pg0220a_raw_0.npz'
rawfile = '/home/srmeeker/scratch/standards/G158-100_raw_0.npz'
#rawfile = '/home/srmeeker/ARCONS-pipeline/examples/Pal2012_landoltPhot/Landolt9542_V_all_0.npz'

#fitfile = '/home/srmeeker/scratch/standards/Landolt9542_fit_0.npz'
#fitfile = '/home/srmeeker/scratch/standards/pg0220a_fit_0.npz'
fitfile = '/home/srmeeker/scratch/standards/G158-100_fit_0.npz'
#fitfile = '/home/srmeeker/ARCONS-pipeline/examples/Pal2012_landoltPhot/Landolt9542_V_fit.npz'


rawDict = np.load(rawfile)
fitDict = np.load(fitfile)

fit = fitDict['fitImg']
params = fitDict['params']
raw = rawDict['stack']
wvls = rawDict['wvls']

print raw

nx,ny = np.shape(raw)[1],np.shape(raw)[2]

xs = np.zeros((nx,ny))
ys = np.zeros((nx,ny))

for i in range(ny):
    xs[:,i] = np.array(np.arange(nx))
for i in range(nx):
    ys[i,:] = np.array(np.arange(ny))

#print xs
#print ys

for n in range(np.shape(raw)[0]):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(xs,ys,fit[n],rstride=1, cstride=1, color='black',alpha=0.1)
    ax.plot_wireframe(xs,ys,raw[n], rstride=1,cstride=1,color='red')

    cset = ax.contour(xs,ys,raw[n], zdir='z', offset=-50, cmap=cm.spring)
    cset = ax.contour(xs,ys,raw[n], zdir='x', offset=-10, cmap=cm.spring)
    cset = ax.contour(xs,ys,raw[n], zdir='y', offset=60, cmap=cm.spring)

    cset = ax.contour(xs,ys,fit[n], zdir='z', offset=-50, cmap=cm.winter)
    cset = ax.contour(xs,ys,fit[n], zdir='x', offset=-10, cmap=cm.winter)
    cset = ax.contour(xs,ys,fit[n], zdir='y', offset=60, cmap=cm.winter)

    ax.set_title(wvls[n])
    ax.set_zlim(0, 170)
    ax.set_xlabel('Col')
    ax.set_ylabel('Row')
    ax.set_zlabel('Counts')

    plt.show()

