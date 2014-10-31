"""
Author: Alex Walter
Date:   Oct 1, 2014

Contains routines for finding hot pixels in a flashing laser cal. 

The supplied ObsFile object (or full filename) is examined. First, 
we find when the light is turned off and on (for each roach) with 
findFlashingTimes(). Then we loop through each pixel and check when 
it goes hot with checkHot(). Finally, an interval time mask is made
indicating when the pixels are hot.

See __main__ for example.

"""
import os, warnings
import numpy as np
from operator import itemgetter
from itertools import groupby
import matplotlib.pyplot as plt
from util.ObsFile import ObsFile
from interval import interval
import headers.TimeMask as tm
import hotpix.hotPixels as hp


def checkHot(avg_on, std_on, avg_off, std_off,onfirst,n_sqrt=2.,n_std=3.,n_photons=1.):
    """
    This code contains the triggers for when a pixel is considered hot. There are several checks:
    1. Check if the count rate error, given by the standard deviation, is consisitent with poisson errors.
    2. Check if neighboring on times have the same count rate (within 'n_std' sigmas of error). Same for off times.
    3. Check if flash on state is brighter than flash off state (at least n_photons brighter).
    4. --NOT IMPLEMENTED YET-- If the pixel is tagged hot for two consecutive off states, make sure the on state between is flagged hot.
    
    Arguments:
    avg_on - list of average count rate during periods of time when the laser is on
    std_on - list of standard deviation 
    avg_off - list of average coutn rate when laser is off
    std_off - list of standard deviations
    onfirst - True if the laser starts out on, False otherwise
    
    Keyword Parameters:
    n_sqrt - Hot if the standard deviation is greater than 'n_sqrt' times the sqrt(average).
    n_std - Hot if the neighboring time periods are more than 'n_std' sigmas different. (Should be a constant count rate)
    n_photons - Hot if on state is less than 'n_photons' brighter than off state.
    
    Careful: The code assumes that the timestream alternates between on and off states. This is important for determining when 
    it goes from on-->off and off-->on. For example, if there are two on states in a row, the comparison to the off state
    will be shifted by one. 
    
    """
    if np.abs(len(avg_on) - len(avg_off)) > 2:
        print "Length of arrays is inconsistent with flashing on and off"
        raise ValueError

    hotmask_on = np.zeros(len(avg_on))
    hotmask_off = np.zeros(len(avg_off))
    
    
    #Check if std is consistent with sqrt
    hotmask_on += std_on > n_sqrt*np.sqrt(avg_on)                   #Should be roughly equal for poison noise
    hotmask_on += std_on == 0                                       #If no standard deviation then flag as hot
    hotmask_off += std_off > n_sqrt*np.sqrt(avg_off)
    std_off[np.where(std_off==0)]=np.sqrt(avg_off)[np.where(std_off==0)]    #Sometimes only 1 bin for off times which makes the std=0
    std_off[np.where(std_off==0)]=1                                         #Sometimes threshold is set high enough that there are no off counts.
    #hotmask_off += std_off == 0
    
    
    #Check if neighboring flashes are consistent
    neighbor_mask = np.abs(avg_on[1:]-avg_on[:-1]) > n_std*np.maximum(std_on[1:],std_on[:-1])       #Should be within error bars
    hotmask_on += np.concatenate(([0],neighbor_mask))                                               #neighbor to right
    hotmask_on += np.concatenate((neighbor_mask,[0]))                                               #neighbor to left
    neighbor_mask = np.abs(avg_off[1:]-avg_off[:-1]) > n_std*np.maximum(std_off[1:],std_off[:-1])
    hotmask_off += np.concatenate(([0],neighbor_mask))
    hotmask_off += np.concatenate((neighbor_mask,[0]))
    

    #Check if flash goes on and off. Have to be careful about if lasers start out being on or off
    length = min(len(avg_on),len(avg_off))                                          #Number of total flashes 
    mask_on2off = (avg_on[:length] - avg_off[:length]) < n_photons                  #Should have at least n_photons more while flashing on
    hotmask_on[:length]+=mask_on2off
    hotmask_off[:length]+=mask_on2off
    #If the laser starts out on then this checks the off-->on transition. Otherwise it checks the on-->off transition. Vice Versa above.
    length = min(len(avg_on[onfirst:]),len(avg_off[not onfirst:]))                  
    mask_off2on = (avg_on[onfirst:length+onfirst] - avg_off[not onfirst:length+(not onfirst)]) < n_photons
    hotmask_on[onfirst:length+onfirst]+=mask_off2on
    hotmask_off[not onfirst:length+(not onfirst)]+=mask_off2on
    
    #Check if time period immediately surrounded by hot periods. 
    #(If the two surrounding off states appear hot then the on state between should be flagged too, & viceversa)
    """
    Add code here
    """
    
    #return boolean arrays. 1 means it's hot during that period of time as determined by the index
    return hotmask_on>0, hotmask_off>0



class flashMask:
    """
    This function creates an exposure time mask for a flashing file. Used with laser cal data where laser is turned on and off

    Provide either the full name of the flashingfile or the ObsFile object
    binSize - determines the timing resolution. (seconds)
              Should divide exposure time evenly
              Should be divided by tickspersec evenly
    flashBuffer - determines a buffer time while flash is transitioning from on-off and off-on. (seconds)
                  Should be an integer number of bins
                  Should be divided by tickspersec evenly
    outputFileName - where to save the h5 file
    """
    def __init__(self, flashingFileName = None, obsFile = None, 
                 binSize = 0.1, flashBuffer = 0.3, 
                 n_sqrt = 2., n_std = 3., n_photons = 1.,startTime=0, endTime= -1,
                 outputFileName = None, verbose=False):

        self.flashingFile=obsFile if obsFile!=None else ObsFile(flashingFileName)
        self.exptime = self.flashingFile.info['exptime']
        self.startTime=startTime if (startTime>=0 and startTime<self.exptime) else 0
        self.endTime=endTime if (endTime>startTime and endTime<=self.exptime) else self.exptime
        
        #Parameters for disecting timestream
        self.flashBuffer_bin = int(round(flashBuffer/binSize))
        self.nBins = int((self.endTime-self.startTime)/binSize)
        self.ticksperbin = int(self.flashingFile.ticksPerSec*binSize)
        
        #Parameters for checking hotness
        self.n_sqrt=n_sqrt
        self.n_std=n_std
        self.n_photons=n_photons
        
        #Parameters for saving
        self.outputFileName=outputFileName
        try:
            os.mkdir(os.path.dirname(outputFileName))
        except:
            pass
        
        #Set verbose to True for output statements
        self.verbose=verbose

    def writeHotPixMasks(self):
        if self.verbose:
            print "Preparing lists of intervals..."
        timeMaskData = []
        for row in range(self.flashingFile.nRow):
            for col in range(self.flashingFile.nCol):
                try:
                    badIntervals = self.convertMasksToBadTimes(row,col)
                    badTimeList = []
                    for eachComponent in badIntervals['laserNotOnIntervals'].components:
                        badTimeList.append((eachComponent[0][0], eachComponent[0][1], tm.timeMaskReason['laser not on'])) #Begin time, end time, flag number  
                    for eachComponent in badIntervals['laserNotOffIntervals'].components:
                        badTimeList.append((eachComponent[0][0], eachComponent[0][1], tm.timeMaskReason['laser not off']))
                    for eachComponent in badIntervals['hotIntervals'].components:
                        badTimeList.append((eachComponent[0][0], eachComponent[0][1], tm.timeMaskReason['hot pixel']))
                    #for eachComponent in badIntervals['unknownIntervals'].components:
                    #    #badTimeList.append((eachComponent[0][0], eachComponent[0][1], tm.timeMaskReason['unknown']))
                    #    badTimeList.append((eachComponent[0][0], eachComponent[0][1], tm.timeMaskReason['hot pixel']))  #Just mask buffer zones as if they were hot
                    badTimeList.sort(key=lambda x: x[0])
                    timeMaskData.append([row, col, badTimeList])
                except IndexError:
                    badTimeList = [(self.startTime*self.flashingFile.ticksPerSec,self.endTime*self.flashingFile.ticksPerSec,tm.timeMaskReason['dead pixel'])]
                    timeMaskData.append([row, col, badTimeList])
        
        if self.verbose:
            print "Writing to "+self.outputFileName
        hp.writeHotPixels(timeMaskData, self.flashingFile, self.outputFileName, startTime=self.startTime, endTime=self.endTime)
        if self.verbose:
            print "\tDone."
        
    def convertMasksToBadTimes(self,row,col):
        masks = self.getPixMasks(row,col)
        laserNotOnMask = -1*(masks['flashMaskOn']*(-1*masks['deadRoachMask']+1))+1      #Mask out these times if we want photons from when laser is on
        laserNotOffMask = -1*(masks['flashMaskOff']*(-1*masks['deadRoachMask']+1))+1    #Mask out these times if we want photons from when laser is off
        hotMask = (masks['hotMaskOn']+masks['hotMaskOff'])>0                            #Mask out these hot times
        
        unknownMask = (masks['flashMaskOn']+masks['flashMaskOff']+masks['deadRoachMask'])==0    #In buffer zone, we don't know if laser is on or off. 
        
        #ticksperbin = int(long(1.0)*self.exptime*self.flashingFile.ticksPerSec/self.nBins)
        ticksperbin = self.ticksperbin
        start=interval([self.startTime*self.flashingFile.ticksPerSec])
        laserNotOnIntervals = [start+interval([i*ticksperbin,(i+1)*ticksperbin]) for i in np.where(laserNotOnMask)[0]]
        laserNotOffIntervals = [start+interval([i*ticksperbin,(i+1)*ticksperbin]) for i in np.where(laserNotOffMask)[0]]
        
        hotIntervals=[start+interval([i*ticksperbin,(i+1)*ticksperbin]) for i in np.where(hotMask)[0]]
        #group successive hot intervals together (along with intervening buffer zone)
        for k, g in groupby(enumerate(np.where((hotMask+unknownMask) > 0)[0]), lambda (i,x):i-x):
            group = map(itemgetter(1), g)
            if len(group) > 2*self.flashBuffer_bin:
                #print group
                hotIntervals.append(start+interval([group[0]*ticksperbin,group[-1]*ticksperbin]))
        #hotIntervals = [start+interval([i*ticksperbin,(i+1)*ticksperbin]) for i in np.where(hotMask)[0]]
        #unknownIntervals = [start+interval([i*ticksperbin,(i+1)*ticksperbin]) for i in np.where(unknownMask)[0]]
        
        #return {'laserNotOnIntervals':interval().union(laserNotOnIntervals), 'laserNotOffIntervals':interval().union(laserNotOffIntervals),
        #        'hotIntervals':interval().union(hotIntervals), 'unknownIntervals':interval().union(unknownIntervals)}
        return {'laserNotOnIntervals':interval().union(laserNotOnIntervals), 'laserNotOffIntervals':interval().union(laserNotOffIntervals),
                'hotIntervals':interval().union(hotIntervals)}
                
        
            
    def getPixMasks(self,row,col):
        r=int(self.flashingFile.beamImage[row,col].split('r')[1][0])
        p=int(self.flashingFile.beamImage[row,col].split('p')[1].split('/')[0])

        k=np.where(np.asarray(self.r_list)==r)[0][0]
        j=np.where(np.asarray(self.pixelNumber_r[k])==p)[0][0]
        
        return {'flashMaskOn':self.flashMaskOn_r[k], 'flashMaskOff':self.flashMaskOff_r[k],
                'hotMaskOn':self.hotMaskOn_r_p[k][j], 'hotMaskOff':self.hotMaskOff_r_p[k][j],
                'deadRoachMask':self.deadRoachMask_r[k]}
        
            
    def findFlashingTimes(self,showPlot=False):
        if self.verbose:
            print "Finding when each roach is flashing..."
        r_list = []
        timestream_r = []
        pixelNumber_r = []
        #Add up timestreams for all pixels in each roach
        for row in range(self.flashingFile.nRow):
            for col in range(self.flashingFile.nCol):
                r=int(self.flashingFile.beamImage[row,col].split('r')[1][0])
                p=int(self.flashingFile.beamImage[row,col].split('p')[1].split('/')[0])
                temp_times = (self.flashingFile.getTimedPacketList(row,col,firstSec=self.startTime, integrationTime= self.endTime-self.startTime))['timestamps']
                temp_timestream,bins = np.histogram(temp_times,self.nBins)

                if len(temp_times)==0:
                    continue

                if (not (r in r_list)):
                    r_list.append(r)
                    timestream_r.append(temp_timestream)
                    pixelNumber_r.append([p])
                else:
                    timestream_r[np.where(np.asarray(r_list)==r)[0]]+=temp_timestream
                    pixelNumber_r[np.where(np.asarray(r_list)==r)[0]].append(p)

        flashMaskOn_r = []
        flashMaskOff_r = []
        deadRoachMask_r = []
        foundFlash=[]
        #Find when lasers are on and off or a timeing error occured
        for k in range(len(r_list)):
            avg = np.average(timestream_r[k])
            mask = timestream_r[k]>avg
            mask_raw = mask
            
            #remove mask times ON that are less than buffer size
            startsOn = np.where((1.0*mask[1:]-mask[:-1])==1)[0]+1
            s_mask = np.array_split(mask,startsOn)
            nOn = np.asarray([np.sum(arr) for arr in s_mask])
            #print k
            #print startsOn
            #print nOn
            s_mask2 = np.asarray(s_mask)
            s_mask2[np.where(nOn<self.flashBuffer_bin)[0]] = np.asarray(s_mask)[np.where(nOn<self.flashBuffer_bin)[0]]*0
            mask = np.concatenate(s_mask2.tolist())
            #remove mask time OFF that are less than buffer size
            startsOff = np.where((1.0*mask[1:]-mask[:-1])==-1)[0]+1
            s_mask = np.array_split(mask,startsOff)
            nOff = np.asarray([np.sum(-1*arr+1) for arr in s_mask])
            s_mask2 = np.asarray(s_mask)
            s_mask2[np.where(nOff<self.flashBuffer_bin)[0]] = np.asarray(s_mask)[np.where(nOff<self.flashBuffer_bin)[0]]*0+1
            mask = np.concatenate(s_mask2.tolist())
            
            flashMaskOn = mask*np.concatenate((mask[self.flashBuffer_bin:],[0]*self.flashBuffer_bin))*np.concatenate(([0]*self.flashBuffer_bin,mask[:-1*self.flashBuffer_bin]))
            mask=-1*mask+1
            flashMaskOff = mask*np.concatenate((mask[self.flashBuffer_bin:],[0]*self.flashBuffer_bin))*np.concatenate(([0]*self.flashBuffer_bin,mask[:-1*self.flashBuffer_bin]))
            
            mask = timestream_r[k]!=0
            deadRoachMask = mask*np.concatenate((mask[self.flashBuffer_bin:],[1]*self.flashBuffer_bin))*np.concatenate(([1]*self.flashBuffer_bin,mask[:-1*self.flashBuffer_bin]))
            deadRoachMask=-1*deadRoachMask+1
            
            flashMaskOn_r.append(flashMaskOn)
            flashMaskOff_r.append(flashMaskOff)
            deadRoachMask_r.append(deadRoachMask)

            if np.sum(flashMaskOn)>len(flashMaskOn)/10. and np.sum(flashMaskOff)>len(flashMaskOff)/100.: 
                foundFlash.append(True)
            else: 
                foundFlash.append(False)
            if showPlot:
                plt.figure()
                #plt.plot(diff,label=r_list[k])
                plt.plot(timestream_r[k],label=str(r_list[k]))
                #plt.plot(filtered_timestream[k],label='filtered')
                #plt.plot(mask*1000-500+avg,label='mask')
                plt.plot((flashMaskOn-0.5)*800+avg,label='on')
                plt.plot((flashMaskOff-0.5)*600+avg,label='off')
                plt.plot((deadRoachMask-0.5)*400+avg,label='dead')
                plt.plot((mask_raw-0.5)*1000+avg,label='raw')
                plt.axhline(y=avg,c='k',ls='--',label='avg')
                plt.legend()
                #plt.show()
        if showPlot:
            plt.show()
        
        self.r_list = r_list
        self.pixelNumber_r = pixelNumber_r
        self.flashMaskOn_r = flashMaskOn_r
        self.flashMaskOff_r = flashMaskOff_r
        self.deadRoachMask_r = deadRoachMask_r
        
        if np.sum(foundFlash)==0.:
            print "Did not observe a flash!"
            raise ValueError
            


    def findHotPixels(self):
        if self.verbose:
            print "Looking for hot pixels..."
    
        hotMaskOn_r_p = []
        hotMaskOff_r_p = []
        
        t = self.flashingFile.beamImage[0,0].split('t')[1]
        for k in range(len(self.r_list)):
            flashMaskOn = self.flashMaskOn_r[k]*(-1*self.deadRoachMask_r[k]+1)
            flashMaskOff = self.flashMaskOff_r[k]*(-1*self.deadRoachMask_r[k]+1)
            
            startsOn = np.where((flashMaskOn[1:]-flashMaskOn[:-1])==1)[0]+1
            #print '\t'+str(len(startsOn))
            #endsOn = np.where((flashMaskOn[1:]-flashMaskOn[:-1])==-1)[0]
            startsOff = np.where((flashMaskOff[1:]-flashMaskOff[:-1])==1)[0]+1
            #endsOff = np.where((flashMaskOff[1:]-flashMaskOff[:-1])==-1)[0]

            hotMaskOn_r_p.append([])
            hotMaskOff_r_p.append([])
            for p in range(len(self.pixelNumber_r[k])):
                pixelLabel = '/r'+str(self.r_list[k])+'/p'+str(self.pixelNumber_r[k][p])+'/t'+str(t)
                [[row],[col]] = np.where(self.flashingFile.beamImage==pixelLabel)
                temp_times = (self.flashingFile.getTimedPacketList(row,col,firstSec=self.startTime, integrationTime= self.endTime-self.startTime))['timestamps']
                temp_timestream,bins = np.histogram(temp_times,self.nBins)

                #Find list of average count rates for each period the laser is on
                split = np.array_split(temp_timestream,startsOn)[1:]  #The first array in split never contains a laser on period. (It's either a time buffer, off period, or empty)
                split_mask = np.array_split(flashMaskOn,startsOn)[1:]
                #split_m = np.ma.array(split,mask=split_mask*-1+1)
                #avg_on = np.ma.getdata(np.ma.mean(split_m,axis=1))
                #std_on = np.ma.getdata(np.ma.std(split_m,axis=1))
                #It'd be nice if we could do this without looping...
                avg_on = np.ma.array([np.ma.mean(np.ma.array(split[i],mask=split_mask[i]*-1+1)) for i in range(len(split))]).filled(0)
                std_on = np.ma.array([np.ma.std(np.ma.array(split[i],mask=split_mask[i]*-1+1)) for i in range(len(split))]).filled(0)

                
                #Find list of average count rates for each period the laser is off
                #split = np.array_split(temp_timestream,startsOff)[1:]  #The first array in split never contains an interesting period of time
                #split_mask = np.array_split(flashMaskOff,startsOff)[1:]
                ##split_m = np.ma.array(split,mask=split_mask*-1+1)  
                ##avg_off = np.ma.getdata(np.ma.mean(split_m,axis=1))
                ##std_off = np.ma.getdata(np.ma.std(split_m,axis=1))
                split = np.array_split(temp_timestream,startsOn)    #By splitting on startsOn trigger we ensure that off states alternate with on states. 
                split_mask = np.array_split(flashMaskOff,startsOn)
                
                if len(split[0])==0 or np.all(split_mask[0]==0):    #A little tricky on the indexing
                    split = split[1:]
                    split_mask = split_mask[1:]
                if np.all(split_mask[-1]==0):
                    split = split[:-1]
                    split_mask = split_mask[:-1]
                #If there are two on states in a row then there will be an avg=0, std=0 filler between for the off state
                avg_off = np.ma.array([np.ma.mean(np.ma.array(split[i],mask=split_mask[i]*-1+1)) for i in range(len(split))]).filled(0)
                std_off = np.ma.array([np.ma.std(np.ma.array(split[i],mask=split_mask[i]*-1+1)) for i in range(len(split))]).filled(0)
                
                
                #####Check if hot#####
                onfirst = startsOn[0]<startsOff[0]
                hotmask_on, hotmask_off = checkHot(avg_on, std_on, avg_off, std_off,onfirst=onfirst, n_sqrt=self.n_sqrt,n_std=self.n_std,n_photons=self.n_photons)

                #Make hot pix masks
                hotpixmask_on = np.array_split(np.zeros(len(flashMaskOn)),startsOn)
                hotpixmask_on = np.transpose(np.add(np.transpose(hotpixmask_on),np.concatenate(([0],hotmask_on)))).tolist()
                hotpixmask_on = np.concatenate(hotpixmask_on)*flashMaskOn
                
                hotpixmask_off = np.array_split(np.zeros(len(flashMaskOff)),startsOn)   #Split as before, on startsOn
                if len(hotpixmask_off[0])==0 or np.all(np.array_split(flashMaskOff,startsOn)[0]==0): hotmask_off = np.concatenate(([0],hotmask_off))   #Fix up indexing mess from before
                if np.all(np.array_split(flashMaskOff,startsOn)[-1]==0): hotmask_off = np.concatenate((hotmask_off,[0]))
                hotpixmask_off = np.transpose(np.add(np.transpose(hotpixmask_off),hotmask_off)).tolist()
                hotpixmask_off = np.concatenate(hotpixmask_off)*flashMaskOff


                #plt.figure
                #plt.plot(temp_timestream,'b-',label='r'+str(self.r_list[k])+' p' +str(p))
                #plt.plot(startsOn,[10]*len(startsOn),'go',label='Starts On')
                #plt.plot(endsOn,[0]*len(endsOn),'go',label='Ends On')
                #plt.plot(startsOff,[10]*len(startsOff),'ro',label='Starts Off')
                #plt.plot(endsOff,[0]*len(endsOff),'ro',label='Ends Off')
                #plt.plot(10*hotpixmask_off,'k-',label='Hot mask')
                #plt.legend()
                #plt.show()


                hotMaskOn_r_p[k].append(hotpixmask_on)
                hotMaskOff_r_p[k].append(hotpixmask_off)
                
        self.hotMaskOn_r_p = hotMaskOn_r_p
        self.hotMaskOff_r_p = hotMaskOff_r_p


if __name__ == '__main__':
    from util.FileName import FileName
    filename = "/ScienceData/PAL2014/20141020/cal_20141021-052251.h5"
    outputFileName = FileName(obsFile=ObsFile(filename)).timeMask()
    #outputFileName = None
    masker = flashMask(filename,startTime = 0,endTime=-1,outputFileName = outputFileName,verbose=True)
    masker.findFlashingTimes(showPlot=False)
    masker.findHotPixels()
    masker.writeHotPixMasks()


