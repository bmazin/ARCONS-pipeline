#This code uses Popup to plot a 3D image

from functools import partial
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from util.popup import *


def plot3DImage(fig,ax,data,fit=None,surface=False,countour=False):
    data[np.where(np.invert(data>0.))]=0.

    x=np.tile(range(len(data[0])),(len(data),1))
    y=np.tile(range(len(data)),(len(data[0]),1)).transpose()
    
    ax = fig.add_subplot(111, projection='3d')
    if surface:
        ax.plot_surface(x, y, data, rstride=1, cstride=1, color='black',alpha=0.1)
    else:
        ax.bar3d(x.flatten(),y.flatten(),x.flatten()*0,1,1,data.flatten(),color='cyan',alpha=0.9)
    
    if countour:
        #cset = ax.contourf(x, y, data, zdir='z', offset=0, cmap=cm.coolwarm)
        cset = ax.contourf(x, y, data, zdir='x', offset=0., cmap=cm.coolwarm)
        cset = ax.contourf(x, y, data, zdir='y', offset=0., cmap=cm.coolwarm)
    else:
        ax.plot(range(len(data[0])+1),np.zeros(len(data[0])+1),zs=np.concatenate(([0],np.amax(data,0))),zdir='z',color='blue',drawstyle='steps')
        ax.plot(np.zeros(len(data)+1),range(len(data)+1),zs=np.concatenate(([0],np.amax(data,1))),zdir='z',color='blue',drawstyle='steps')

    if fit!=None:
        ax.plot_wireframe(x+0.5, y+0.5, fit, rstride=1, cstride=1, color='red')
        ax.plot(np.asarray(range(len(data[0])))+0.5,np.zeros(len(data[0])),zs=np.amax(fit,0),zdir='z',color='red')
        ax.plot(np.zeros(len(data)),np.asarray(range(len(data)))+0.5,zs=np.amax(fit,1),zdir='z',color='red')

    ax.set_zlim(0, np.amax(data))
    ax.set_xlim(0,np.amax(x))
    ax.set_ylim(0,np.amax(y))
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
