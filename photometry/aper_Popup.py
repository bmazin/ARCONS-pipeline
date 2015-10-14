
#LightCurve Popup for checking lightcurve, aperture photometry, interpolations...
from functools import partial
import numpy as np
import matplotlib
from util import utils
from util.popup import *
from photometry.LightCurve import LightCurve

def clickCanvas(self,LC,event):
    if event.inaxes is self.axes and self.mpl_toolbar._active is None:
        time = LC.photometry_dict['startTimes']
        closest_time_ind = np.argmin(np.abs(time - event.xdata))
        #print event.xdata, ' --> ', time[closest_time_ind]
        
        pop=PopUp(parent=self,title='JD: '+str(time[closest_time_ind])+' 2D Image')
        pop_2DImage(pop,LC,closest_time_ind)
        pop.show()
        
        pop=PopUp(parent=self,title='JD: '+str(time[closest_time_ind])+' Interpolated Image')
        pop_interpolatedImage(pop,LC,closest_time_ind)
        pop.show()
        
        print '\nimage:',closest_time_ind
        #print 'centroid:',LC.centroids[closest_time_ind]
        


        
def pop_interpolatedImage(self,LC,ind):
    interpolation = LC.photometry_dict['interpolation'][ind]
    aper_radius = LC.photometry_dict['apertureRad'][ind]
    annulus_inner = LC.photometry_dict['annulusInnerRad'][ind]
    annulus_outer = LC.photometry_dict['annulusOuterRad'][ind]
    
    image = LC.im_dict['images'][ind]
    centroids = LC.centroids[ind]
    interpImage = utils.interpolateImage(image, method=interpolation)
    angleArray = np.linspace(0,2*np.pi,180)
    
    self.plotArray(image=interpImage, title='Interpolated Image',cmap=matplotlib.cm.gnuplot2)
    for star_i in range(len(centroids)):
        self.axes.plot(centroids[star_i][0]+aper_radius[star_i]*np.cos(angleArray), centroids[star_i][1] + aper_radius[star_i]*np.sin(angleArray), 'b-')
        self.axes.plot(centroids[star_i][0]+annulus_inner[star_i]*np.cos(angleArray), centroids[star_i][1] + annulus_inner[star_i]*np.sin(angleArray), 'r-')
        self.axes.plot(centroids[star_i][0]+annulus_outer[star_i]*np.cos(angleArray), centroids[star_i][1] + annulus_outer[star_i]*np.sin(angleArray), 'r-')
        
        self.axes.plot(centroids[star_i][0],centroids[star_i][1],'gx')
    
    
        
def pop_2DImage(self,LC,ind):

    interpolation = LC.photometry_dict['interpolation'][ind]
    aper_radius = LC.photometry_dict['apertureRad'][ind]
    annulus_inner = LC.photometry_dict['annulusInnerRad'][ind]
    annulus_outer = LC.photometry_dict['annulusOuterRad'][ind]
    
    image = LC.im_dict['images'][ind]
    image[np.invert(np.isfinite(image))]=0.
    centroids = LC.centroids[ind]
    angleArray = np.linspace(0,2*np.pi,180)
    
    self.plotArray(image=image, title='Image',cmap=matplotlib.cm.gnuplot2)
    for star_i in range(len(centroids)):
        self.axes.plot(centroids[star_i][0]+aper_radius[star_i]*np.cos(angleArray), centroids[star_i][1] + aper_radius[star_i]*np.sin(angleArray), 'b-')
        self.axes.plot(centroids[star_i][0]+annulus_inner[star_i]*np.cos(angleArray), centroids[star_i][1] + annulus_inner[star_i]*np.sin(angleArray), 'r-')
        self.axes.plot(centroids[star_i][0]+annulus_outer[star_i]*np.cos(angleArray), centroids[star_i][1] + annulus_outer[star_i]*np.sin(angleArray), 'r-')
        
        self.axes.plot(centroids[star_i][0],centroids[star_i][1],'gx')


def hoverCanvas(self,time,event):
    if event.inaxes is self.axes:
        closest_time_ind = np.argmin(np.abs(time - event.xdata))
        self.status_text.setText(str(time[closest_time_ind]))

def plotLightCurve(self,LC):
    LC.loadLightCurve(photometryType='aper')
    time = LC.photometry_dict['startTimes']
    flags = LC.photometry_dict['flag']
    flux=LC.photometry_dict['flux']
    
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
    
def plotLightCurveRatio(self,LC):
    LC.loadLightCurve(photometryType='aper')
    time = LC.photometry_dict['startTimes']
    flags = LC.photometry_dict['flag']
    flux=LC.photometry_dict['flux']
    
    targetFlux=flux[:,0]
    for star in range(len(flux[0])-1):
        lbl=LC.targetName + '/Ref'+str(star)
        ref_flux = flux[:,star+1]
        self.axes.plot(time[np.where(flags==0.)],targetFlux[np.where(flags==0.)]/ref_flux[np.where(flags==0.)],'.',label=lbl)

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
    identifier = '15s_4000-9000A_flat_hp_V3'
    
    LC=LightCurve(fileID=identifier,path=path,targetName=None,run=None,verbose=True,showPlot=False) 
    
    pop(plotFunc = partial(plotLightCurve,LC=LC),title="Aper Light Curve Popup")
    pop(plotFunc = partial(plotLightCurveRatio,LC=LC),title="Aper Light Curve Popup")










