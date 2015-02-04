
#LightCurve Popup for checking lightcurve, aperture photometry, interpolations...
from functools import partial
import numpy as np
import matplotlib
from util import utils
from util.popup import *
from photometry.LightCurve import LightCurve

def clickCanvas(self,LC,event):
    if event.inaxes is self.axes and self.mpl_toolbar._active is None:
        photometryDict = LC.getLightCurve(photometryType='aperture',fromPhotometryFile=True,save=False)
        time = photometryDict['startTimes']
        closest_time_ind = np.argmin(np.abs(time - event.xdata))
        #print event.xdata, ' --> ', time[closest_time_ind]
        
        pop=PopUp(parent=self,title='JD: '+str(time[closest_time_ind])+' 2D Image')
        pop_2DImage(pop,LC,closest_time_ind)
        pop.show()
        
        pop=PopUp(parent=self,title='JD: '+str(time[closest_time_ind])+' Interpolated Image')
        pop_interpolatedImage(pop,LC,closest_time_ind)
        pop.show()


        
def pop_interpolatedImage(self,LC,ind):
    #These should really be saved in the photometry file and grabbed from there...
    interpolation = "linear"
    aper_radius = 5
    annulus_inner = 10
    annulus_outer = 15
    
    image = LC.im_dict['images'][ind]
    centroid = LC.centroids[ind]
    interpImage = utils.interpolateImage(image, method=interpolation)
    angleArray = np.linspace(0,2*np.pi,180)
    
    self.plotArray(image=interpImage, title='Interpolated Image',cmap=matplotlib.cm.gnuplot2)
    for star_i in range(len(centroid)):
        self.axes.plot(centroid[star_i][0]+aper_radius*np.cos(angleArray), centroid[star_i][1] + aper_radius*np.sin(angleArray), 'b-')
        self.axes.plot(centroid[star_i][0]+annulus_inner*np.cos(angleArray), centroid[star_i][1] + annulus_inner*np.sin(angleArray), 'r-')
        self.axes.plot(centroid[star_i][0]+annulus_outer*np.cos(angleArray), centroid[star_i][1] + annulus_outer*np.sin(angleArray), 'r-')
        
    
    
        
def pop_2DImage(self,LC,ind):
    image = LC.im_dict['images'][ind]
    image[np.invert(np.isfinite(image))]=0.

    self.plotArray(image=image, title='Image',cmap=matplotlib.cm.gnuplot2)

def hoverCanvas(self,time,event):
    if event.inaxes is self.axes:
        closest_time_ind = np.argmin(np.abs(time - event.xdata))
        self.status_text.setText(str(time[closest_time_ind]))

def plotLightCurve(self,LC):
    photometryDict = LC.getLightCurve(photometryType='aperture',fromPhotometryFile=True,save=False)
    time = photometryDict['startTimes']
    flags = photometryDict['flag']
    flux=photometryDict['flux']
    
    for star in range(len(flux[0])):
        lbl=LC.targetName
        if star>0: lbl='ref'+str(star-1)
        star_flux = flux[:,star]
        self.axes.plot(time[np.where(flags==0.)],star_flux[np.where(flags==0.)],'.',label=lbl)
        #self.axes.plot(time[np.where(flags==1.)],star_flux[np.where(flags==1.)],'r.')
        #self.axes.plot(time[np.where(flags>1.1)],star_flux[np.where(flags>1.1)],'ro')
    #self.axes.plot(time[np.where(flags==1.)],star_flux[np.where(flags==1.)],'r.',label='Centroid Fail')
    #self.axes.plot(time[np.where(flags>1.1)],star_flux[np.where(flags>1.1)],'ro',label='Fit Fail')
    self.axes.legend()
    self.axes.set_xlabel('Julian Date')
    self.axes.set_ylabel('Integrated Stellar Flux [counts/sec]')
    
    cid = self.fig.canvas.mpl_connect('button_press_event',partial(clickCanvas,self,LC))
    cid = self.fig.canvas.mpl_connect('motion_notify_event', partial(hoverCanvas,self,time))
    
    self.fig.canvas.draw()

if __name__ == '__main__':
    #path = '/Scratch/DisplayStack/PAL2014/HAT_P1'
    #identifier = '1'
    path = '/Scratch/DisplayStack/PAL2014/1SWASP_J2210'
    identifier = '0'
    
    LC=LightCurve(fileID=identifier,path=path,targetName=None,run=None,verbose=True,showPlot=False) 
    
    pop(plotFunc = partial(plotLightCurve,LC=LC),title="Light Curve Popup")










