#Flashing Wavecal viewer
import numpy as np
from util.popup import *
from util.ObsFile import ObsFile
from util.FileName import FileName
from functools import partial
from flashMask import *
import smooth

def clickCanvas(laserFlashFile,self,event):
    

    if event.inaxes is self.axes:
        col = int(round(event.xdata))
        row = int(round(event.ydata))
        print "Clicked: ("+str(col)+', '+str(row)+') --> '+laserFlashFile.beamImage[row,col]
        
        if not hasattr(self,'laserFlashFile'):
            self.laserFlashFile = laserFlashFile
            #print 'Making masks...'
            #self.masker = flashMask(obsFile=laserFlashFile)
            #self.masker.findFlashingTimes()
            #self.masker.findHotPixels()
            #print '\tDone.'
        
        pop=PopUp(parent=self,title='Timestream')
        pop.laserFlashFile = laserFlashFile
        pop_timestream(pop,row,col)
        pop.show()
        
        pop=PopUp(parent=self,title='Spectra')
        pop.laserFlashFile=laserFlashFile
        pop_spectra(pop,row,col)
        pop.show()

def pop_spectra(self,row,col):
    self.laserFlashFile.switchOffHotPixTimeMask()
    dataDict = self.laserFlashFile.getTimedPacketList(row,col,timeSpacingCut=0.001)
    peakHeights=np.asarray(dataDict['peakHeights'])*1.0
    baselines=np.asarray(dataDict['baselines'])*1.0
    peakHeights-=baselines
    biggest_photon = int(min(peakHeights))
    n_inbin,phase_bins=np.histogram(peakHeights,bins=np.abs(biggest_photon),range=(biggest_photon,0))
    phase_bins=(phase_bins+(phase_bins[1]-phase_bins[0])/2.0)[:-1]
    try:
        last_ind = np.where(n_inbin>5)[0][-1]
    except IndexError:
        last_ind=len(n_inbin)-1
    expTime = self.laserFlashFile.hotPixTimeMask.expTime
    self.axes.plot(phase_bins,n_inbin*1.0/expTime, 'k.',alpha=0.5,label="raw")
    self.axes.set_xlim(phase_bins[(np.where(n_inbin >= 3))[0][0]],phase_bins[last_ind])

    self.laserFlashFile.switchOnHotPixTimeMask(reasons=['laser not on','hot pixel'])
    dataDict = self.laserFlashFile.getTimedPacketList(row,col,timeSpacingCut=0.001)
    peakHeights=np.asarray(dataDict['peakHeights'])*1.0
    baselines=np.asarray(dataDict['baselines'])*1.0
    peakHeights-=baselines
    biggest_photon = int(min(peakHeights))
    n_inbin,phase_bins=np.histogram(peakHeights,bins=np.abs(biggest_photon),range=(biggest_photon,0))
    phase_bins_1=(phase_bins+(phase_bins[1]-phase_bins[0])/2.0)[:-1]
    intTime_1 = self.laserFlashFile.hotPixTimeMask.getEffIntTime(row,col)
    smoothed_1 = np.asarray(smooth.smooth(n_inbin, 30, 'hanning'))
    self.axes.plot(phase_bins_1,n_inbin*1.0/intTime_1, 'b.',label="laser on, hotpix masked")
    
    self.laserFlashFile.switchOnHotPixTimeMask(reasons=['laser not off','hot pixel'])
    dataDict = self.laserFlashFile.getTimedPacketList(row,col,timeSpacingCut=0.001)
    peakHeights=np.asarray(dataDict['peakHeights'])*1.0
    baselines=np.asarray(dataDict['baselines'])*1.0
    peakHeights-=baselines
    biggest_photon = int(min(peakHeights))
    n_inbin,phase_bins=np.histogram(peakHeights,bins=np.abs(biggest_photon),range=(biggest_photon,0))
    phase_bins_2=(phase_bins+(phase_bins[1]-phase_bins[0])/2.0)[:-1]
    intTime_2 = self.laserFlashFile.hotPixTimeMask.getEffIntTime(row,col)
    smoothed_2 = np.asarray(smooth.smooth(n_inbin, 30, 'hanning'))
    self.axes.plot(phase_bins_2,n_inbin*1.0/intTime_2, 'g.',label="laser off, hotpix masked")
    
    self.axes.legend(loc=2)
    
    self.axes.plot(phase_bins_1, smoothed_1*1.0/intTime_1,'b-')
    self.axes.plot(phase_bins_2,smoothed_2*1.0/intTime_2,'g-')
    
    self.axes.set_xlabel('phase [ADC/DAC units]')
    self.axes.set_ylabel('Exposure time adjusted count rate [photons/s]')


def pop_timestream(self,row,col):
    self.laserFlashFile.switchOffHotPixTimeMask()
    times = self.laserFlashFile.getTimedPacketList(row,col)['timestamps']
    factor = 10
    exptime = self.laserFlashFile.info['exptime']
    nBins=exptime*factor
    timestream,bins = np.histogram(times,bins=nBins,range=(0,exptime))
    self.axes.plot(timestream,'k-',alpha=0.5,label='raw')

    self.laserFlashFile.switchOnHotPixTimeMask(reasons=['laser not on'])
    times_on = self.laserFlashFile.getTimedPacketList(row,col)['timestamps']
    #nBins=3000
    timestream_on,bins = np.histogram(times_on,bins=nBins,range=(0,exptime))
    self.axes.plot(timestream_on,'b-',label='laser on')

    self.laserFlashFile.switchOnHotPixTimeMask(reasons=['laser not off'])
    times_off = self.laserFlashFile.getTimedPacketList(row,col)['timestamps']
    #nBins=3000
    timestream_off,bins = np.histogram(times_off,bins=nBins,range=(0,exptime))
    self.axes.plot(timestream_off,'g-',label='laser off')
    
    bad_time_intervals=self.laserFlashFile.getPixelBadTimes(row,col,reasons=['hot pixel'])
    for eachComponent in bad_time_intervals.components:
        #print eachComponent[0]
        self.axes.axvspan(eachComponent[0][0]*factor,eachComponent[0][1]*factor,facecolor='r',alpha=0.5)
        
    self.axes.axvspan(-2,-1,facecolor='r',alpha=0.5,label='hot')
    #print eachComponent[0][0],eachComponent[0][1]
    self.axes.legend()

    self.axes.set_xlabel('time [s/'+str(factor)+']')
    self.axes.set_ylabel('count rate [photons/s]')


if __name__ == '__main__':
    filename = "/ScienceData/PAL2014/20141020/cal_20141021-052251.h5"
    laserFlashFile = ObsFile(filename)
    print 'Loading hotpix file'
    laserFlashFile.loadHotPixCalFile(FileName(obsFile=laserFlashFile).timeMask(), reasons=['hot pixel','laser not on'])
    print 'Making image'
    #result = laserFlashFile.getPixelCountImage(integrationTime= 2,getRawCount=True)
    #image = result['image']
    #image = result['effIntTimes']
    image = laserFlashFile.hotPixTimeMask.getEffIntTimeImage()
    print 'Making plot'
    plotArray(image,button_press_event=partial(clickCanvas,laserFlashFile))




