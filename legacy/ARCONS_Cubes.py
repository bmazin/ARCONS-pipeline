'''
Created on Aug 12, 2011

@author: kobrien
Aug 12, 2011: added h5cubes class to interact with cubes (x,y,Energy) 
              of histogrammed data
Dec 14, 2011: added h5lights class to interact with lightcurves
Jun 11, 2012: fixed bug in peakfit
'''
from tables import *
import numpy as np
import struct
import time
#from decimal import *
import os
import subprocess
from subprocess import Popen, PIPE, STDOUT
import ephem, datetime
import pyfits

def __init__(self):
    pass

class h5histtable(IsDescription):
    """The pytables derived class that hold histogram'ed data on the disk
    """
    # recorded data
    hist = Float32Col(2000)

def peakfit(y1,y2,y3):
    y4=y2-0.125*((y3-y1)**2)/(y3+y1-2*y2)
    return y4

def peakfittime(y1,y2,y3):
    y4=np.float128(0.125)*(y1-y3)/(y3+y1-np.float128(2.0)*y2)
    return y4

class h5lights():
    def genlight(self,filename,mask,tstart=-1,tend=-1,peakfitbool=True):
        '''
        generates lightcurve based on masked regions between a start and stop time
        It then sorts the resulting photon time stream.
        It will fit the peak of the 3 samples stored in the packet and use the 
        height and time of the resulting fit, otherwise it uses the middle point.
        It returns the time and pulseheight.  
        self.jd0 is jd at start of lightcurve
        v1.0 Kieran O'Brien - Dec 2011
        '''
        h5file_temp=openFile(filename)
        self.exptime=int(h5file_temp.root.header.header.col('exptime')[0])
        if(tstart >= self.exptime):
            print ' ERROR: start time after end of file'
            return
        if(tend > self.exptime):
            print ' ERROR: end time after end of file, leave blank for tend=exptime'
            return
        if(tstart < 0):
            tstart=0
        if(tend < 0):
            tend=int(h5file_temp.root.header.header.col('exptime')[0])
        if(tend-tstart == 0):
            print ' ERROR:  no times to load in range:', tstart, tend
            return
#        try:
        print ' Generating lightcurve...'
#        print h5file_temp.root.header.header.col('jd')[0]
#        print h5file_temp.root.header.header.col('ut')[0]
#        self.mjdrefi = int(h5file_temp.root.header.header.col('ut')[0]/86400.0)
#        self.mjdreff = (h5file_temp.root.header.header.col('ut')[0]-86400.0*self.mjdrefi)/86400.0
#        self.mjdrefi = 40587+self.mjdrefi
#        print self.mjdrefi, self.mjdreff
#        print np.float128(self.mjdreff.astype('float128')/86400.0)
#        print np.float128((0.000001+self.mjdreff)/86400.0)
#        print h5file_temp.root.header.header.col('localtime')[0]
        self.jd0=h5file_temp.root.header.header.col('jd')[0].astype('float128')
        self.mjdrefi=np.floor(self.jd0-np.float128(2400000.5))
        self.mjdreff=np.float128(self.jd0)-np.float128(2400000.5)-np.float128(self.mjdrefi)
        print self.mjdrefi,self.mjdreff
        masked=np.where(mask > 0.5)
        bmap=h5file_temp.root.beammap.beamimage.read()
        obstimes=np.array([],dtype=np.float128)
        obsheights=np.array([],dtype=np.float128)
        for i in range(np.shape(masked)[1]):
            for k in range(tstart,tend):
                photons= h5file_temp.root._f_getChild(bmap[masked][i])[k]
                if(peakfitbool):
                    obsheights_y1= np.right_shift(photons,20)%4096
                    obsheights_y2= np.right_shift(photons,32)%4096
                    obsheights_y3= np.right_shift(photons,44)%4096
                    obsheights_y1=np.array(obsheights_y1,dtype=np.float128)
                    obsheights_y2=np.array(obsheights_y2,dtype=np.float128)
                    obsheights_y3=np.array(obsheights_y3,dtype=np.float128)
                    obstimes1 = np.float128(k)+np.float128(1.0e-6)*(np.array(photons%1048576,dtype=np.float128)+np.float128(peakfittime(obsheights_y1,obsheights_y2,obsheights_y3)))
                    obsheights1 = peakfit(obsheights_y1,obsheights_y2,obsheights_y3)
#                    print np.array(photons%1048576,dtype=np.float128)+np.float128(0.000001)*np.float128(peakfittime(obsheights_y1,obsheights_y2,obsheights_y3))
                else:
                    obstimes1=np.float128(k)+np.float128(photons%1048576)/np.float128(1.0e+6)
                    obsheights1= np.right_shift(photons,32)%4096
                obstimes=np.append(obstimes,obstimes1)
#            print obstimes1[0],obstimes[0]
                obsheights=np.append(obsheights,obsheights1)
        # sort the array in time
        idx=np.argsort(obstimes)
        self.obstimes=obstimes[idx]
        self.obsheights=obsheights[idx]
        print ' ... finished'
#        except:
#            print " ERROR: problem generating lightcurve"
        h5file_temp.close()
        return
    
    def savelight(self,filename,bary=False):
        '''
        outputs a file containing the times and pulse heights
        v1.0 Kieran O'Brien - Dec 2011
        '''
        try:
            outfile=open(filename,'w')
            if bary:
                print 'saving barytime corrected lightcurve'
                for i in range(len(self.barytimes)):
                    outfile.write(str("%.17f" %self.barytimes[i])+'\t'+str("%.17f" %self.obsheights[i])+'\n')
            else:
                print 'saving uncorrected lightcurve'
                for i in range(len(self.obstimes)):
                    outfile.write(str("%.17f" %self.obstimes[i])+'\t'+str("%.17f" %self.obsheights[i])+'\n')
            outfile.close()
        except:
            print 'ERROR: problem writing file'
        return

    def savefits(self,filename,bary=False):
        '''
        outputs a FITS file containing the histogram lightcurve
        v1.0 Kieran O'Brien - Dec 2011
        '''
#        try:
        if bary:
            print 'saving barytime corrected lightcurve'
            col1=pyfits.Column(name='BARYTIME',format='E', array=self.barytimes)
            col2=pyfits.Column(name='COUNTS',format='E', array=self.obsheights)
            cols=pyfits.ColDefs([col1,col2])
            tbhdu=pyfits.new_table(cols)
            hdu=pyfits.PrimaryHDU(0)
            thdulist=pyfits.HDUList([hdu,tbhdu])
            thdulist.writeto(filename)
        else:
            print 'saving uncorrected lightcurve'
            col1=pyfits.Column(name='TIME',format='E', array=self.obstimes)
            col2=pyfits.Column(name='COUNTS',format='E', array=self.obsheights)
            cols=pyfits.ColDefs([col1,col2])
            tbhdu=pyfits.new_table(cols)
            hdu=pyfits.PrimaryHDU(0)
            thdulist=pyfits.HDUList([hdu,tbhdu])
            thdulist.writeto(filename)
#        except:
#            print 'ERROR: problem writing file'
        return

    def genmask_square(self,xpixel_min,xpixel_max,ypixel_min,ypixel_max):
        '''
        generates a mask for genlight
        v1.0 Kieran O'Brien - Dec 2011
        '''
        self.mask=np.zeros((32,32))
        for i in range(xpixel_min,xpixel_max+1):
            for j in range(ypixel_min,ypixel_max+1):
                self.mask[i,j]=1
        return
    
    def genmask_circ(self,cen,radius):
        '''
        generates a mask for genlight
        v1.0 Kieran O'Brien - Dec 2011
        '''
        self.mask=np.zeros((32,32))
        ivals=np.linspace(0,31,32)
        jvals=np.linspace(0,31,32)
#        print cen[0], cen[1]
        for i in range(32):
            for j in range(32):
                dist=(ivals[i]-cen[0])*(ivals[i]-cen[0])+(jvals[j]-cen[1])*(jvals[j]-cen[1])
#                print i,j,ivals[i],jvals[j],dist,radius*radius
                if(dist < radius*radius):
                    self.mask[i,j]=1
        return

    def foldlight(self,period,nbins,ephemeris=0):
        '''
        folds a lightcurve on a known linear ephemeris
        v1.0 Kieran O'Brien - Dec 2011
        '''
        phase=(self.obstimes-np.float128(ephemeris))/np.float128(period)
        phase=phase-phase.astype('int')
        self.phase,self.phasebins=np.histogram(phase,nbins)
        return
        
    def genbarytimes(self,parfile):
        '''
        Interface to TEMPO2 to create the barycentre corrected lightcurves
        Very slow, but only needs to be run once
        v1.0 Kieran O'Brien - Dec 2011
        '''
        print 'First photon: ', self.obstimes[0]
        print 'Last photon: ',self.obstimes[-1]
        print 'number of photons: ',len(self.obstimes)
        mjdstart = self.mjdrefi + self.mjdreff + self.obstimes[0]/np.float128(86400.0)
        mjdend = self.mjdrefi + self.mjdreff + self.obstimes[-1]/np.float128(86400.0)
        self.barytimes=np.zeros_like(self.obstimes)
        maxlines=np.float128(9999)
        nfiles=1+int(np.floor(len(self.obstimes)/maxlines))
        print nfiles
        for j in range(nfiles):
            jstart=int(j*(maxlines-1))
            jend=int(np.min([len(self.obstimes),jstart+maxlines-1]))
            print jstart+1, ' - ', jend, ' of ',len(self.obstimes)
            # set-up for entry into tempo
            outfile=open('zzz_test.tim','w')
            outfile.write('FORMAT 1\n')
            for i in range(jstart,jend):
                time=np.float128(self.mjdreff)+np.float128(self.obstimes[i])/np.float128(86400.0)
                outfile.write(' w040206_070831.FT 500000000.000 '+str("%5i" % self.mjdrefi)+str(time)[1:21]+' 0.10000 h\n')
            outfile.close()
            # run in tempo and catch the output
            os.putenv('TEMPO2','/usr/share/tempo2')
            proc1=subprocess.Popen(['/usr/share/tempo2/bin/tempo2','-output','general2','-s','"{BAT}"', '-f',parfile,'zzz_test.tim'],stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            proc1_out,proc1_err = proc1.communicate()
            proc1.wait()            
            proc1out=proc1_out.split('"')
            proc1out=np.array(proc1out)
            newtime=np.cast[np.float128](proc1out[(proc1out != '')][1:-1])
            self.barytimes[jstart:jend]=np.float128(86400.)*(newtime-mjdstart)
        return
        
        
class h5cubes():
    def gencube(self,filename):
        '''
        outputs a array containing the histogrammed data (x,y,pulseheight)
        needs self.nhistbins to be set in advance
        v1.0 Kieran O'Brien - Aug 2011
        '''
#        fitpeak=False
        fitpeak=True
        h5file=openFile(filename)
        self.xpix=np.array(32)
        self.ypix=np.array(32)
        try:
            self.bmap = h5file.root.beammap.beamimage.read()
            self.header = h5file.root.header.header.read()
            # make a cube with 32x32x(number of bins in histogram)
            self.tstart=0
            self.tend=int(h5file.root.header.header.col('exptime')[0])

            self.obshist=np.zeros((32,32,self.nhistbins))
#            print 'loading photon data...'
            for i in range(self.xpix):
                for j in range(self.ypix):
#                    print i,j, self.bmap[i][j], self.tstart, self.tend
                    photons= np.concatenate(h5file.root._f_getChild(self.bmap[i][j])[self.tstart:self.tend])
#                    print 'loaded photons...'
                    if(fitpeak):
                        obsheights= peakfit(np.right_shift(photons,20)%4096,np.right_shift(photons,32)%4096,np.right_shift(photons,44)%4096)
                    else:
#                        obsheights= np.right_shift(photons,20)%4096
                        obsheights= np.right_shift(photons,32)%4096
#                        obsheights= np.right_shift(photons,44)%4096
                    self.obshist[i][j],self.histbins = np.histogram(obsheights,bins=self.nhistbins,range=(self.histmin,self.histmin+self.nhistbins))
            self.histbins=0.5*(self.histbins[1:]+self.histbins[:-1])
        except:
            print 'h5 file not complete!'
            pass
        h5file.close()

    def gencube_energy(self,filename,efilename,tstart=0,tend=0):
        '''
        outputs an array containing the histogrammed data (x,y,Energy)
        uses the energyscale files create by Seth
        Stores various calibration files into the file for later use
        v1.0 Kieran O'Brien - Aug 2011
        '''
        fitpeak=True
#        fitpeak=False
        h5file=openFile(filename)
        escalefile=openFile(efilename)
        self.tstart=tstart
        if(tend == 0):
            self.tend=int(h5file.root.header.header.col('exptime')[0])
        else:
            self.tend=np.min([tend,int(h5file.root.header.header.col('exptime')[0])])
        if(self.tstart > self.tend):
            print 'start time after end time'
        self.xpix=np.array(32)
        self.ypix=np.array(32)
        # checking that beam maps match
        bmap1=h5file.root.beammap.beamimage.read()
        bmap2=escalefile.root.beammap.beamimage.read()
        self.bmap = h5file.root.beammap.beamimage.read()
        for i in range(self.xpix):
            for j in range(self.ypix):
                tree1 = bmap1[i][j].split('/')
                tree2 = bmap2[i][j].split('/')
                if (tree1[0] != tree2[0]):
                    if (tree1[1] != tree2[1]):
                        print "Beam maps don\'t match!!!"
                        break
                self.bmap[i][j]=str('/'+tree1[1]+'/'+tree1[2])
    #            print bmap1[i][j],self.bmap[i][j]
                        
        try:
            self.header = h5file.root.header.header.read()
            # make a cube with 32x32x(number of bins in histogram)
            self.escale=escalefile.root.calparams.energyscale
    
            self.obshist=np.zeros((32,32,self.nhistbins))
            for i in range(self.xpix):
                for j in range(self.ypix):
                    photons= np.concatenate(h5file.root._f_getChild(bmap1[i][j])[self.tstart:self.tend])
                    if(fitpeak):
                        obsheights= self.escale[i][j][1]+self.escale[i][j][0]*peakfit(np.right_shift(photons,20)%4096,np.right_shift(photons,32)%4096,np.right_shift(photons,44)%4096)
                    else:
                        obsheights= self.escale[i][j][1]+self.escale[i][j][0]*np.right_shift(photons,32)%4096
#                        print i,j,self.bmap[i][j],self.escale[i][j]
                    self.obshist[i][j],self.histbins = np.histogram(obsheights,bins=self.nhistbins,range=(self.histmin,self.histmax))
            self.histbins=0.5*(self.histbins[1:]+self.histbins[:-1])
        except:
            print 'h5 file not complete!'
            pass
        h5file.close()
        escalefile.close()

    def savehistcube(self, filename,wmode):
        '''
        outputs a file containing the histogtrammed data
        v1.0 Kieran O'Brien - Aug 2011
        '''

        # open the original file and update it
        newh5file = openFile(filename,mode=wmode)

        for i in range(self.xpix):
            for j in range(self.ypix):
                tree = self.bmap[i][j].split('/')
#                print tree
                # if there is no existing roach group, create one
                try:
                    group = newh5file.getNode('/',tree[1] )
                except:
#                    print 'creating roach group'
                    group = newh5file.createGroup("/",tree[1], 'ROACH Identifier' )

                # make a group for each resonator
                try:
                    group = newh5file.getNode(group,tree[2] )
                except:
#                    print 'creating pixel group'
                    group = newh5file.createGroup(group,tree[2],'pixel number' )        

                # put data into table
                try:
                    table = newh5file.createArray(group,'histogram', self.obshist[i][j]) 
                except:
                    print 'failed to create table'
                    pass

                try:
                    table.flush()
                except:
                    print 'failed to flush table'       
        # save the bin values in a separate place
        try:
            table = newh5file.createArray('/','bins', self.histbins) 
            table.flush()
        except:
            print 'failed to create bin position table'
            pass                
        try:
            newh5file.close()
            print 'file closed'
            print 'written: ', filename
        except:
            print 'failed to close file'

    def loadcube(self, filename):
        '''
        Loads a file containing the histogrammed data
        Works on pulse height or energy data
        v1.0 Kieran O'Brien - Aug 2011
        '''
        h5file=openFile(filename,mode='r')
        self.bmap=h5file.root.beammap.beamimage.read()
        self.xpix=np.array(32)
        self.ypix=np.array(32)
        self.exptime=h5file.root.header.header.col('exptime')[0]
        self.histbins=h5file.root.bins.read()
        self.nhistbins=len(self.histbins)
#        print self.nhistbins
        self.obshist=np.zeros((32,32,self.nhistbins))
#        bmapbit=self.bmap.split('/')
        for i in range(32):
            for j in range(32):
                bmapbit=self.bmap[j][i].split('/')
                roachno=bmapbit[1]
                pixelno=bmapbit[2]
                pixelpos=h5file.getNode('/'+roachno, pixelno)
                self.obshist[j][i]=pixelpos.histogram.read()
        try:
            self.qualityflag=h5file.root.calparams.qualityflag.read()
        except:
            print 'failed to load goodpix array'
        try:
            self.ecut=h5file.root.calparams.energycut.read()
        except:
            print 'failed to energy cut-off array'
        h5file.close()


    def rebincube(self,factor):
        '''
        rebins a cube in pulse-height axis
        v1.0 Kieran O'Brien - Aug 2011
        '''

        # rebins cubes, ignores partially filled bins (assumed to be beyond region of interest)
        orighist=self.obshist
        origbins=np.array(self.histbins,dtype=float)
        nbins=int(len(origbins)/factor)
        self.obshist=np.zeros((32,32,nbins))
        self.histbins=np.arange(nbins,dtype=float)
        self.nhistbins=nbins
        for i in range(32):
            for j in range(32):
                for k in range(nbins):
                    self.obshist[i][j][k]=np.sum(orighist[i][j][k*factor:(k+1)*factor])
        for k in range(nbins):
            self.histbins[k]=float(np.mean(origbins[k*factor:(k+1)*factor]))

    def gengoodpix(self):
        '''
        generates a 32x32 array with 1 if good, 0 if bad (eg. no energy scale)
        v1.0 Kieran O'Brien - Aug 2011
        '''
        self.qualityflag=np.zeros((32,32))
        for i in range(32):
            for j in range(32):
#                print np.sum(self.obshist[i][j][:])
                self.qualityflag[i][j]=np.sum(self.obshist[i][j][:])
                self.qualityflag[np.where(self.goodpix > 0.)] = 1
#                print self.goodpix[i][j]
 
    def savegoodpix(self, filename):
        '''
        outputs a file containing the new pixel array
        v1.0 Kieran O'Brien - Aug 2011
        '''
        # open the original file and update it
        newh5file = openFile(filename,mode='a')
        # save the bin values in a separate place
        try:
            newh5file.removeNode('/calparams','qualityflag',recursive=True)
            print 'removed original quality flag map'
        except:
            pass
                
        try:
            filt1 = Filters(complevel=1, complib='blosc', fletcher32=False)
            table = newh5file.createCArray('/calparams', 'qualityflag',Int32Atom(), (32,32), filters=filt1)
            table[:,:] = self.qualityflag[:,:]
#            table = newh5file.createArray('/calparams','qualityflag', self.qualityflag) 
            newh5file.flush()
        except:
            print 'quality flag map not saved'
            pass
                
        try:
            newh5file.close()
            print 'file closed'
            print 'written: ', filename
        except:
            print 'failed to close file'

    def noisemodel(self):
        '''
        Generates an array based on a simple photon noise model
        v1.0 Kieran O'Brien - Aug 2011
        '''
        # noise model has sqrt of number of photons, apart from where histogram is zero
        self.noise=np.sqrt(self.obshist)
        # find low energy cut-off
        cutoff=np.zeros((self.obshist.shape[0],self.obshist.shape[1]),dtype=int)
        for i in range(self.obshist.shape[0]):
            for j in range(self.obshist.shape[1]):
#                self.noise[i][j]=np.max([1.,self.noise[i][j]])
#                print self.qualityflag[i][j]
#                print self.obshist[i][j][:]
                if(self.qualityflag[i][j] > 0. and np.sum(self.obshist[i][j][:]) > 0.):
                    noflux=np.where(self.obshist[i][j] == 0.)[0]
                    noflux=np.append(noflux,[self.nhistbins+2])
                    if(len(noflux) > 1 ):
#                        print noflux, noflux[1:], noflux[:-1]
#                    print np.where(noflux[0][1:-1]-noflux[0][:-2] >1.)[0]
                        gaps=np.where(noflux[1:]-noflux[:-1] >1.)[0]
#                        print gaps
                        cutoff[i][j]=gaps[0]
                        cutoff[i][j]=np.min([cutoff[i][j]+5,self.nhistbins])
                    else:
                        cutoff[i][j]=self.nhistbins
                        print 'here', noflux
#                    print noflux
#                    print i,j,cutoff[i][j]
                else:
                    cutoff[i][j]=self.nhistbins
                # increase the noise below the cut-off energy (no data)
                self.noise[i][j][:cutoff[i][j]]=1e+9     
                
    def energy2idx(self, energy):
        try:
            idx=np.where(self.histbins < energy)[0][-1]
        except:
            idx=0
        return idx
        
def slice_image(input,start,end): 
        '''
        makes an image from a region of the cube
        v1.0 Kieran O'Brien - Aug 2011
        '''
        slice=np.zeros((input.shape[0],input.shape[1]))
        for i in range(input.shape[0]):
            for j in range(input.shape[1]):
                slice[i,j]=np.sum(input[i,j,start:end])
        return slice
       
                    
def swap_beammap(infile, newbeammap, outfile, photon=False): 
        '''
        does what is says on the tin
        v1.0 Kieran O'Brien - Aug 2011
        '''
        inh5file=openFile(infile,mode='r')
        inh5file.copyFile(outfile,overwrite=True)
        bmap=openFile(newbeammap,mode='r')
        outh5file=openFile(outfile,mode='a')
        outh5file.copyNode(bmap.root.beammap,newparent=outh5file.root,recursive=True,overwrite=True)         
        # catch if the file is a photon file and needs the extension to the beammap
        if(photon):
            beamimage=np.asarray(bmap.root.beammap.beamimage)
            timestamp=inh5file.root.r0.p0._v_children.keys()[0]
            for i in range(32):
                for j in range(32):
                    beamimage[i][j]=str(bmap.root.beammap.beamimage[i][j]+timestamp)
            outh5file.removeNode('/beammap','beamimage')
            filt1 = Filters(complevel=1, complib='blosc', fletcher32=False)
            table = outh5file.createCArray('/beammap', 'beamimage',StringAtom(itemsize=40), (32,32), filters=filt1)
            table[:,:] = beamimage[:,:]
            outh5file.flush()
        inh5file.close()
        outh5file.close()
        bmap.close()

def copyHeader(headfile,outfile,wmode):
        '''
        outputs a file containing only the header, useful as a starting place 
        for new files
        v1.0 Kieran O'Brien - Aug 2011
        '''

        infile=openFile(headfile,mode='r')
        outfile=openFile(outfile,mode=wmode)
        try:
            outfile.copyNode(infile.root.header,newparent=outfile.root,recursive=True)  
        except:
            print 'no header'
            pass
        outfile.close()        
        infile.close()

def updateHeaderTime(infile,outfile):
        '''
        outputs a file containing the corrected times based on the start time 
        of the observation and not the time the file was written
        v1.0 Kieran O'Brien - Aug 2011
        '''
        inh5file=openFile(infile,mode='r')
        inh5file.copyFile(outfile,overwrite=True)
        inh5file.close()
        outh5file=openFile(outfile,mode='a')
#        try:
        timestamp=int(outh5file.root.r0.p0._v_children.keys()[0].strip('t'))
        dt=datetime.datetime.utcfromtimestamp(timestamp)
        jd=np.float64(ephem.julian_date(dt))
        table=outh5file.root.header.header
        table.cols.jd[0]=jd
        table.cols.ut[0]=timestamp
        table.flush()    
#            outfile.copyNode(infile.root.header,newparent=outfile.root,recursive=True)  
#        except:
#            print 'no header'
#            pass
        outh5file.close()        

def copyBeammap(beamfile,outfile,wmode):
        infile=openFile(beamfile,mode='r')
        outfile=openFile(outfile,mode=wmode)
        try:
            outfile.copyNode(infile.root.beammap,newparent=outfile.root,recursive=True)         
        except:
            print 'no beammap'
            pass
        outfile.close()        
        infile.close()

def copyEscale(escalefile,outfile,wmode):
        infile=openFile(escalefile,mode='r')
        outfile=openFile(outfile,mode=wmode)
        try:
            outfile.copyNode(infile.root.calparams,newparent=outfile.root,recursive=True)         
        except:
            print 'no calib info copied'
            pass
        outfile.close()        
        infile.close()
        
def ev2nm(energy):
        return 1240./energy

def nm2ev(wave):
        return 1240./wave

def repeatphase(phibin,phi):
        return np.append(phibin,1.0+phibin),np.append(phi,phi)
    
def foldlight(light,period,nbins,ephemeris=0,pdot=0):
        '''
        folds array on a given ephemeris.
        v1.0 Kieran O'Brien - Jun 2012
        '''
        if(pdot == 0):
            phi=(light-np.float128(ephemeris))/np.float128(period)
            phi=phi-phi.astype('int')
            phase,phasebins=np.histogram(phi,nbins)
        else:
            deltat=(light-np.float128(ephemeris))
            phi=deltat/(period+pdot*deltat)
            phi=phi-phi.astype('int')
            phase,phasebins=np.histogram(phi,nbins)
        return phasebins, phase
'''
def foldlight(light,period,nbins,ephemeris=0,pdot=0,variance=True):
        if(pdot == 0):
            if(variance):
                phi=(light-np.float128(ephemeris))/np.float128(period)
                phi=phi-phi.astype('int')
                for i in range(nbins):
                    idx=np.where((phi >= phimin)&&(phi < phimax))[0]
                    phase[i]=np.mean(phi[idx])
            
            else:
                phi=(light-np.float128(ephemeris))/np.float128(period)
                phi=phi-phi.astype('int')
                phase,phasebins=np.histogram(phi,nbins)
        else:
            deltat=(light-np.float128(ephemeris))
            phi=deltat/(period+pdot*deltat)
            phi=phi-phi.astype('int')
            phase,phasebins=np.histogram(phi,nbins)
        return phasebins, phase
'''
def foldlight_freq(light,freq,nbins,ephemeris=0,freqdot=0):
        if(freqdot == 0):
            phi=(light-np.float128(ephemeris))*np.float128(freq)
            phi=phi-phi.astype('int')
            phase,phasebins=np.histogram(phi,nbins)
        else:
            deltat=(light-np.float128(ephemeris))
            phi=deltat*(freq+freqdot*deltat)
            phi=phi-phi.astype('int')
            phase,phasebins=np.histogram(phi,nbins)
        return phasebins, phase
        