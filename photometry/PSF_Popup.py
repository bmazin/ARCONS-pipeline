
#LightCurve Popup for checking lightcurve, PSF fits, etc. 

from functools import partial
import numpy as np
import matplotlib
from util.popup import *
from util.fitFunctions import model_list
from photometry.LightCurve import LightCurve
from photometry.plot3DImage import plot3DImage


def clickCanvas(self,LC,event):
    if event.inaxes is self.axes and self.mpl_toolbar._active is None:
        photometryDict = LC.getLightCurve(fromLightCurveFile=True)
        time = photometryDict['startTimes']
        closest_time_ind = np.argmin(np.abs(time - event.xdata))
        print event.xdata, ' --> ', time[closest_time_ind]
        
        pop=PopUp(parent=self,title='JD: '+str(time[closest_time_ind])+' PSF fit')
        pop_image(pop,LC,closest_time_ind)
        pop.show()
        
        pop=PopUp(parent=self,title='JD: '+str(time[closest_time_ind])+' Fit Residual')
        pop_residualImage(pop,LC,closest_time_ind)
        pop.show()
        
        pop=PopUp(parent=self,title='JD: '+str(time[closest_time_ind])+' 2D Image')
        pop_2DImage(pop,LC,closest_time_ind)
        pop.show()
        
def pop_image(self,LC,ind,model='multiple_2d_circ_gauss_func'):
    photometryDict = LC.getLightCurve(fromLightCurveFile=True)
    im_dict = LC.getImages(fromImageStack=True)
    
    image = im_dict['images'][ind]
    image[np.invert(np.isfinite(image))]=0.
    errs = np.sqrt(image)
    errs[np.where(image==0.)]=np.inf
    parameters = photometryDict['parameters'][ind]
    models = model_list[model](parameters)(p=np.ones(len(parameters)),data=image,return_models=True)
    guess = models[0]
    for m in models[1:]:
        guess+=m
    
    plot3DImage(self.fig,self.axes,image,errs=errs,fit=guess)
    
def pop_residualImage(self,LC,ind,model='multiple_2d_circ_gauss_func'):
    photometryDict = LC.getLightCurve(fromLightCurveFile=True)
    im_dict = LC.getImages(fromImageStack=True)
    
    image = im_dict['images'][ind]
    image[np.invert(np.isfinite(image))]=0.
    errs = np.sqrt(image)
    errs[np.where(image==0.)]=np.inf
    parameters = photometryDict['parameters'][ind]
    models = model_list[model](parameters)(p=np.ones(len(parameters)),data=image,return_models=True)
    guess = models[0]
    for m in models[1:]:
        guess+=m
    
    residualImage = (image-guess)/np.sqrt(image)
    residualImage[np.where(image==0)] = 0
    residualImage[np.invert(np.isfinite(residualImage))]=0.
    self.plotArray(image=residualImage, title='Image Residual',cmap=matplotlib.cm.gnuplot2)
    
def pop_2DImage(self,LC,ind):
    im_dict = LC.getImages(fromImageStack=True)
    
    image = im_dict['images'][ind]
    image[np.invert(np.isfinite(image))]=0.

    self.plotArray(image=image, title='Image',cmap=matplotlib.cm.gnuplot2)
    
def hoverCanvas(self,time,event):
    if event.inaxes is self.axes:
        closest_time_ind = np.argmin(np.abs(time - event.xdata))
        self.status_text.setText(str(time[closest_time_ind]))

def plotLightCurve(self,LC):
    photometryDict = LC.getLightCurve(fromLightCurveFile=True)
    time = photometryDict['startTimes']
    flags = photometryDict['flag']
    flux=photometryDict['flux']
    
    for star in range(len(flux[0])):
        lbl='target'
        if star>0: lbl='ref'+str(star-1)
        star_flux = flux[:,star]
        self.axes.plot(time[np.where(flags==0.)],star_flux[np.where(flags==0.)],'.-',label=lbl)
        self.axes.plot(time[np.where(flags==1.)],star_flux[np.where(flags==1.)],'r.')
        self.axes.plot(time[np.where(flags>1.1)],star_flux[np.where(flags>1.1)],'ro')
    self.axes.plot(time[np.where(flags==1.)],star_flux[np.where(flags==1.)],'r.',label='Centroid Fail')
    self.axes.plot(time[np.where(flags>1.1)],star_flux[np.where(flags>1.1)],'ro',label='Fit Fail')
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
    LC = LightCurve(path,fileID = identifier, PSF=True) 
    
    pop(plotFunc = partial(plotLightCurve,LC=LC),title="Light Curve Popup")









