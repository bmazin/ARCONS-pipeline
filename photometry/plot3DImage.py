#This code uses Popup to plot a 3D image

from functools import partial
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from util.popup import *


def plot3DImage(fig,ax,data,fit=None):
    data[np.where(np.invert(data>0.))]=0.

    x=np.tile(range(len(data[0])),(len(data),1))
    y=np.tile(range(len(data)),(len(data[0]),1)).transpose()
    
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(x, y, data, rstride=1, cstride=1, color='black',alpha=0.1)
    #cset = ax.contourf(x, y, data, zdir='z', offset=0, cmap=cm.coolwarm)
    cset = ax.contourf(x, y, data, zdir='x', offset=0, cmap=cm.coolwarm)
    cset = ax.contourf(x, y, data, zdir='y', offset=np.amax(y), cmap=cm.coolwarm)
    if fit!=None:
        ax.plot_wireframe(x, y, fit, rstride=1, cstride=1, color='red')

    ax.set_zlim(0, np.amax(data))
    cid = fig.canvas.mpl_connect('scroll_event', partial(scroll3D,fig,ax))
    
def scroll3D(fig,ax,event):
    increment = 0.05
    currentZlim = ax.get_zlim()[1]
    if event.button == 'up':
        newZlim = currentZlim-increment*currentZlim
    if event.button == 'down':
        newZlim = currentZlim+increment*currentZlim
        
    if newZlim < 10:
        newZlim=10
        
    ax.set_zlim(0,newZlim)
    fig.canvas.draw()
