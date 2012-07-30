#!/usr/bin/env python
# encoding: utf-8
"""
CalFit.py

Fits 5-peak pixels very well, often still fails on 4-peaks.
Updated from CalFit_092111 to not include 0 count bins in fit.

Created by Seth Meeker on 2011-09-06.
Copyright (c) 2011 . All rights reserved.
"""

import numpy as np
import time
from tables import *
import matplotlib
#matplotlib.use('PDF')
import matplotlib.pylab as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
import mpfit

def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

def BIRFUNC(x,p):
# fit parameters for wavelength cal with blue+IR laser
# 0 = height of primary blue peak
# 1 = location of primary blue peak
# 2 = sigma of blue primary peak
# 3 = linear term of dispersion
# 4 = constant term of dispersion
# 5 = scale factor for sigma
# 6 = scale factor for Blue peaks
# 7 = height of IR peak

    peaksev = [3.102,2.516,1.925,1.273]

    #position of B1 peak
    B1pos = p[1]
    #sigma of B1 peak
    B1sig = p[2]
    #height of B1 peak
    B1height = p[0]
    
    z1 = (x-B1pos)/B1sig
    bp1 = B1height*np.exp((-z1**2)/2.0)

    #position of B2 peak
    B2pos = (peaksev[1]-p[4])/p[3]
    #sigma of B2 peak
    B2sig = p[5]*B1sig
    #height of B2 peak
    B2height = p[6]*B1height
    
    z2 = (x-B2pos)/B2sig
    bp2 = B2height*np.exp((-z2**2)/2.0)

    #position of B3 peak
    B3pos = (peaksev[2]-p[4])/p[3]
    #sigma of B3 peak
    B3sig = p[5]*B1sig
    #height of B3 peak
    B3height = p[6]*B2height
    
    z3= (x-B3pos)/B3sig   
    bp3 = B3height*np.exp((-z3**2)/2.0)

    #position of IR peak
    IRpos = (peaksev[3]-p[4])/p[3]
    #sigma of IR peak
    IRsig = B1sig
    #height of IR peak
    IRheight = p[7]
    
    z4 = (x-IRpos)/IRsig  
    irp1 = IRheight*np.exp((-z4**2)/2.0)

    fit = bp1 + bp2 + bp3 + irp1
    
    return fit 

def GAUSSFUNC(x,p):
# fit parameters for wavelength cal with blue+IR laser
# 0 = height of primary blue peak
# 1 = location of primary blue peak
# 2 = sigma of blue primary peak

    #position of B1 peak
    B1pos = p[1]
    #sigma of B1 peak
    B1sig = p[2]
    #height of B1 peak
    B1height = p[0]
    
    z1 = (x-B1pos)/B1sig
    bp1 = B1height*np.exp((-z1**2)/2.0)
    
    return bp1

def BIRCALFIT(p, fjac=None, x=None, y=None, err=None):
    fit = BIRFUNC(x,p)    
    status=0
    #print len(y)
    #print len(fit)
    #print len(err)
    return [status, (y-fit)/err]

def GAUSSFIT(p,fjac=None, x=None, y=None, err=None):
    fit = GAUSSFUNC(x,p)    
    status=0    
    return [status, (y-fit)/err]

def onclick(event):
    global flag
    print "Stored Good Pixel"
    flag = 1



#############################  BEGIN MAIN PROGRAM  ##############################


if __name__ == '__main__':
    
    peaksev = [3.102,2.516,1.925,1.273]
    cutoffev = 1.12
    
    #pp = PdfPages('./223423/cal_223423_full.pdf')
    #matplotlib.rcParams['font.size']=6
    
    # Load up configuration file
    configfile = './calibration.cfg'
    
    outfile = "./045533/CalFitOut_all_full.txt"
    pfile = "./045533/FitParams_all_full.txt"
    #outfile = "./223423/junk_full.txt"
    #pfile = "./223423/junk_params_full.txt"
    
    obslist = []
    pixlist = []
    skiplist = []
    skipped=0
    
    cfg = open(configfile,'r')
    for line in cfg:
        label,junk,data = line.partition('=')
        if label == 'pixel':
            #pass
            pixel=data.strip()
            #pixlist.append(pixel)
            #print pixel
        elif label == 'dir':
            dir=data.strip()
            print 'Working in directory ' + dir
        else:
            obslist.append(label.strip())
            #file = label.strip()
            print file
    
    file = obslist[0]
    
    h5file = openFile(file, mode = "r")  
    ts = int(h5file.root.header.header.col('ut')[0] + 2.0)
    exptime = int(h5file.root.header.header.col('exptime')[0])
    target = h5file.root.header.header.col('target')[0]
    
    
    #to use list of pixels in config file, comment out this block and uncomment the append in the cfg file read above
    bmap = h5file.root.beammap.beamimage.read()
    for i in range(32):
        for j in range(32):
            pstring = bmap[i][j]
            pixel, timejunk = pstring.split('t')
            pixel.strip()
            #if int(pixel[2]) == 2 or int(pixel[2]) == 3:
            pixlist.append(pixel)
    
    f=open(outfile, 'w')
    f.close()
    
    fp = open(pfile, 'w')
    fp.write("pixel\tflag\tp0low\tp0val\tp0high\tp1low\tp1val\tp1high\tp2low\tp2val\tp2high\tp3low\tp3val\tp3high\tp4low\tp4val\tp4high\tp5low\tp5val\tp5high\tp6low\tp6val\tp6high\tp7low\tp7val\tp7high\n")
    fp.close()
    
    for pixel in pixlist:
        f=open(outfile, 'a')
        fp = open(pfile, 'a')
        
        pn = pixel+'t%d' % (ts)
        data= np.concatenate(h5file.root._f_getChild(pn)[:])
        
        print '\n----------------------------------------------\n', pn
        
        if len(data)==0:
            print "Dead pixel"
            
            continue
        
        caldata = np.right_shift(data,32)%4096
        counts, bins = np.histogram(caldata, 400, range = (0,2000))
        
        print (bins[1]-bins[0])/2.0
        bins += (bins[1]-bins[0])/2.0
        
        #smooth histogram and locate first blue peak
        wlen=9
        smoothed = smooth(counts,window_len=wlen,window='flat')
        
        fig=plt.figure(1)
        ax=fig.add_subplot(111)
        ax.plot(bins[:-1],counts)
        ax.set_title(pn)
        
        global flag
        flag = 0
        
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        
        #ax.plot(bins[:-1],smoothed)
        
        derivs = [0]
        for i in range(len(smoothed)-1):
            derivs.append((smoothed[i+1]-smoothed[i])/(bins[i+1]-bins[i]))

        derivs = np.array(derivs)
        #print derivs
        
        #get indices for derivs
        #startderivs = (abs(derivs) > np.mean(abs(derivs))).nonzero()[0]
        startderivs = (abs(derivs) > 0.0).nonzero()[0]
        #check that peak location found is above 50 counts; need some limit to ensure it doesnt find a noise peak
        startcounts = (counts > 50).nonzero()[0]
        
        #find where in the startderivs list is also in startcounts
        for m in range(len(startderivs)):
            if startderivs[m] in startcounts:
                start = startderivs[m]
                break
        
        newderivs = derivs[start:]
        
        try:
            index = (newderivs < 0).nonzero()[0][0]
        except IndexError:
            skipped+=1
            skiplist.append(pixel)
            print "1 Skipped pixel ", pixel
            m=0
            b=0
            flag=0
            f.write("%9s \t %7f \t %6f \t%1d \n" %(pn,m,b,flag))
            f.close()
            fp.close()
            ax.clear()
            continue
        
        B1ind = index+start-1
        
        B1peak = bins[B1ind]
        B1height = smoothed[B1ind]
        
        #print B1ind
        #print derivs[B1ind]
        #print bins[B1ind]
        #print smoothed[B1ind]
        #print counts[B1ind]
        
        ax.plot(np.repeat(bins[B1ind],B1height),np.arange(int(B1height)))
        
        '''
        fig2=plt.figure(2)
        ax2=fig2.add_subplot(111)
        ax2.plot(bins[:-1], derivs)
        plt.show()
        '''
        
        #IR peak height appears loosely related to B1 peak height
        IRheight = B1height*0.7753+291
        
        tail = counts[:B1ind]
        hmind = (tail > (B1height/2.0)).nonzero()[0][0]
        tailind = (tail > (B1height/4.0)).nonzero()[0][0]
        #print (B1height/2.0)
        #print 'hmind =', hmind
        hmpos = bins[hmind]
        sigbins = B1ind-hmind
        #print 'B1bin =',B1peak
        #print 'hmpos =',hmpos
        B1sig = (B1peak-hmpos)
        
        
        ######------------
        #counts = smoothed
        ######------------
        
        print "First B1sig = ", B1sig
        
        ###First initial fit to get better estimate of sigma of first peak
        
        # 0 = height of primary blue peak
        # 1 = location of primary blue peak
        # 2 = sigma of blue primary peak
        parinfo=[ {'value':0., 'fixed':0, 'limited':[1,1], 'limits':[0.,0.]} for i in range(3) ]
        
        parinfo[0]['value'] = B1height #+ B1height*np.random.normal(scale=0.2) 
        parinfo[0]['limits'] = [counts[B1ind-sigbins], B1height+100]
        parinfo[0]['fixed'] = 0
        
        parinfo[1]['value'] = B1peak #+ B1peak*np.random.normal(scale=0.1) 
        parinfo[1]['limits'] = [B1peak-(B1sig/2.0), B1peak+(B1sig/2.0)]
        parinfo[1]['fixed'] = 0
    
        parinfo[2]['value'] = B1sig #+ (B1sig)*np.random.normal(scale=0.1) 
        parinfo[2]['limits'] = [B1sig/2.0,B1sig*1.5]
        parinfo[2]['fixed'] = 0
        
        g1y = counts[:B1ind]
        g2y = np.flipud(g1y)
        gy = np.append(g1y,g2y)
        gy = gy[:(len(bins)-1)]
        gx = bins[:-1]
        sig0 = np.sqrt(gy)
        sig0[:tailind]=1.0e9
        sig0[(sigbins+B1ind):]=1.0e9

        fa = {'x':gx, 'y':gy, 'err':sig0}  
        try:
            m = mpfit.mpfit(GAUSSFIT,functkw=fa,parinfo=parinfo,quiet=1)
            p0 = m.params
            B1sig = p0[2]
        #print m.status,m.fnorm
        except ValueError:
            skipped+=1
            skiplist.append(pixel)
            print "Error estimating B1sig for", pixel
            #B1sig = 
            #ax.clear()
            #m=0
            #b=0
            #flag=0
            #f.write("%9s \t %7f \t %6f \t%1d \n" %(pn,m,b,flag))
            #f.close()
            #fp.close()
            #continue
        
        #m appears to correlate nicely with B1sig
        
        #new m and b relations
        minverse = -3.8729 * B1sig - 50.869
        mguess = 1.0/minverse
        bguess = peaksev[0]-mguess*B1peak
        
        #mguess = 0.00008513*B1sig - 0.00839664
        #bguess = peaksev[0]-mguess*B1peak
        
        print "New B1sig = ", B1sig
        print "mguess = ", mguess
        
        
        '''
        fig2 = plt.figure(2)
        ax2=fig2.add_subplot(111)
        ax2.plot(gx,gy)
        fit0 = GAUSSFUNC(gx,p0)
        ax2.plot(gx,fit0)
        ax2.plot(np.repeat(gx[hmind],1000),np.arange(1000))
        ax2.plot(np.repeat(gx[sigbins+B1ind],1000),np.arange(1000))
        '''
        
        ax.plot(np.repeat(hmpos,int(counts[hmind])),np.arange(int(counts[hmind])))
        #plt.show()
        
        #only give fit nonzero data
        try:
            first0=(counts[:B1ind] <= 0).nonzero()[0][-1]
        except IndexError:
            first0=0
        print 'index of blue 0 data cutoff = ',first0
        
        try:
            last0=(counts[B1ind:] <=0).nonzero()[0][0] 
            last0 += B1ind
        except IndexError:
            last0 = len(counts)-1
        print 'index of IR 0 data cutoff = ',last0
        
        fitx =bins[first0:last0]
        fity =counts[first0:last0]
        
        #ax.plot(fitx,fity)
        #plt.show()
        
        #set up mpfit params and fit
        bestfnorm = 1e20
        
        print "fitting pixel ", pn
        for i in xrange(1):
            
            # 0 = height of primary blue peak
            # 1 = location of primary blue peak
            # 2 = sigma of blue primary peak
            # 3 = linear term of dispersion
            # 4 = constant term of dispersion
            # 5 = scale factor for sigma
            # 6 = scale factor for Blue peaks
            # 7 = height of IR peak
            parinfo=[ {'value':0., 'fixed':0, 'limited':[1,1], 'limits':[0.,0.]} for i in range(8) ]
            
            parinfo[0]['value'] = B1height #+ B1height*np.random.normal(scale=0.2) 
            parinfo[0]['limits'] = [counts[B1ind-sigbins], B1height+100]
            parinfo[0]['fixed'] = 0
            
            parinfo[1]['value'] = B1peak #+ B1peak*np.random.normal(scale=0.1) 
            parinfo[1]['limits'] = [B1peak-(B1sig/2.0), B1peak+(B1sig/2.0)]
            parinfo[1]['fixed'] = 0
        
            parinfo[2]['value'] = B1sig #+ (B1sig)*np.random.normal(scale=0.1) 
            parinfo[2]['limits'] = [B1sig-B1sig*0.2,B1sig+B1sig*0.2]
            parinfo[2]['fixed'] = 0
        
            parinfo[3]['value'] = mguess #+ mguess*np.random.normal(scale=0.1)  
            parinfo[3]['limits'] = [mguess-0.0015,mguess+0.0015]
            #parinfo[3]['limits'] = [-0.0075,-0.003]
            parinfo[3]['fixed'] = 0
        
            parinfo[4]['value'] = bguess #+ bguess*np.random.normal(scale=0.15)
            parinfo[4]['limits'] = [bguess-3,bguess+3]       
            
            parinfo[5]['value'] = 1.1 #+ 1*np.random.normal(scale=0.1)
            parinfo[5]['limits'] = [0.5,1.5]
            
            parinfo[6]['value'] = (0.80) #+ 0.66*np.random.normal(scale=0.4)
            parinfo[6]['limits'] = [0.66,1.2]    
        
            parinfo[7]['value'] = IRheight #+ B1height*np.random.normal(scale=0.2)
            parinfo[7]['limits'] = [0.5*IRheight,1.3*IRheight]

            sigma = np.sqrt(fity)/2.0
            for k in range(len(sigma)):
                if sigma[k] == 0:
                    sigma[k]=1
                    #sigma[k]=1.0e9
                if fity[k] > counts[B1ind]:
                    sigma[k]=1.0e9
                    
            #do not try to fit tail to the left of half max of blue peak
            #sigma[:tailind] = 1.0e9
            
            #print len(fitx)
            #print len(fity)
            #print len(sigma)
            fa = {'x':fitx, 'y':fity, 'err':sigma}
            #fa = {'x':bins[:-1], 'y':counts, 'err':sigma}    
            
            #print parinfo
            
            #cutoffpos = (cutoffev-parinfo[4]['value'])/parinfo[3]['value']
            #cutoffind = (bins > cutoffpos).nonzero()[0][0] - 1 
            #ax.plot(np.repeat(cutoffpos,1000),np.arange(1000))
            #sigma[cutoffind:]=1.0e9
            
            print "mpfitting"
            m = mpfit.mpfit(BIRCALFIT,functkw=fa,parinfo=parinfo,quiet=1)
            p = m.params
            print parinfo
            print m.status,m.fnorm
            if( m.status > 0 and m.fnorm < bestfnorm):
                bestfnorm = m.fnorm
                bestp = p
            
        p1 = bestp  
        print "First fit p= ",
        print p1
        
        fit1 = BIRFUNC(bins[:-1],p1)
        ax.plot(bins[:-1],fit1)
        #plt.show()
        
        '''
        #estimate IR peak position to set cutoff threshold for 5th peak
        IRfitpeak = (peaksev[3]-p1[4])/p1[3]
        IRfitind = (bins > IRfitpeak).nonzero()[0][0]
        IRhmfitind = (fit1[IRfitind:] < fit1[IRfitind]*0.75 ).nonzero()[0][0]
        IRhmfitind+=IRfitind
        print "IRfitind=", IRfitind
        print "IRhmfitind=", IRhmfitind 
        ax.plot(np.repeat(bins[IRhmfitind],counts[IRhmfitind]),np.arange(counts[IRhmfitind]))
        sigma[IRhmfitind:]=1.0e9
        '''
        
        #sigma = np.sqrt(fity)/2.0
        #for k in range(len(sigma)):
           #if sigma[k] == 0:
                #sigma[k]=1.0
                #sigma[k]=1.0e9
                
        #do not try to fit tail to the left of half max of blue peak
        #sigma[:tailind] = 1.0e9
        
        #instead of IR peak, use first fit and a cutoff energy to set cutoff threshold for 5th peak
        cutoffpos = (cutoffev-p1[4])/p1[3]
        try:
            cutoffind = (fitx > cutoffpos).nonzero()[0][0] - 1
        except IndexError:
            cutoffind = len(fitx)-2
        
        ax.plot(np.repeat(cutoffpos,300),np.arange(300))
        sigma[cutoffind:]=1.0e9
        
        #for k in range(len(sigma[IRhmfitind:])):
                #if counts[k] == 0:
                    #sigma[k]=1.0
                    #sigma[k]=1.0e9
        
        for i in xrange(1):
            
            # 0 = height of primary blue peak
            # 1 = location of primary blue peak
            # 2 = sigma of blue primary peak
            # 3 = linear term of dispersion
            # 4 = constant term of dispersion
            # 5 = scale factor for sigma
            # 6 = scale factor for Blue peaks
            # 7 = height of IR peak
            parinfo=[ {'value':0., 'fixed':0, 'limited':[1,1], 'limits':[0.,0.]} for i in range(8) ]
            
            parinfo[0]['value'] = B1height #+ B1height*np.random.normal(scale=0.2) 
            parinfo[0]['limits'] = [counts[B1ind-sigbins], B1height+100]
            parinfo[0]['fixed'] = 0
            
            parinfo[1]['value'] = B1peak #+ B1peak*np.random.normal(scale=0.05) 
            parinfo[1]['limits'] = [B1peak-(B1sig/2.0), B1peak+(B1sig/2.0)]
            parinfo[1]['fixed'] = 0
        
            parinfo[2]['value'] = p[2] #+ B1sig*np.random.normal(scale=0.05) 
            parinfo[2]['limits'] = [p[2]-p[2]*0.1,p[2]+p[2]*0.1]
            parinfo[2]['fixed'] = 0
        
            parinfo[3]['value'] = (p1[3]) #+ (p1[3])*(np.random.normal(scale=0.05))    
            parinfo[3]['limits'] = [p1[3]-0.003,p1[3]+0.003]
            #parinfo[3]['limits'] = [-0.007,-0.003]
        
            parinfo[4]['value'] = (p1[4]) #+ (p1[4])*(np.random.normal(scale=0.05))
            parinfo[4]['limits'] = [p1[4]*0.5,p1[4]*1.5]       
            
            parinfo[5]['value'] = (p1[5]) #+ (p1[5])*(np.random.normal(scale=0.05))
            parinfo[5]['limits'] = [p1[5]*0.66,p1[5]*1.33]
            
            parinfo[6]['value'] = p1[6] #+ p1[6]*np.random.normal(scale=0.1)
            parinfo[6]['limits'] = [0.66,1.2]    
        
            parinfo[7]['value'] = p1[7] #+ p1[7]*np.random.normal(scale=0.2)
            parinfo[7]['limits'] = [0.8*p1[7],1.3*p1[7]]  

            #sigma = np.sqrt(counts)
            #for k in range(len(sigma)):
                #if sigma[k] == 0:
                   #sigma[k]=1.0
            #do not try to fit tail to the left of half max of blue peak
            #sigma[:tailind] = 1.0e9
            
            fa = {'x':fitx, 'y':fity, 'err':sigma}
            
            #print parinfo
            
            m = mpfit.mpfit(BIRCALFIT,functkw=fa,parinfo=parinfo,quiet=1)
            p = m.params
            print m.status,m.fnorm
            if( m.status > 0 and m.fnorm < bestfnorm):
                bestfnorm = m.fnorm
                bestp = p
        
        p2 = bestp  
        print "Second fit p= ",
        print p2
        
        fit2 = BIRFUNC(bins[:-1],p2)
        ax.plot(bins[:-1],fit2)

        cutoffpos = (cutoffev-p2[4])/p2[3]
        try:
            cutoffind = (fitx > cutoffpos).nonzero()[0][0] - 1
        except IndexError:
            cutoffind = len(fitx)-2
        
        ax.plot(np.repeat(cutoffpos,300),np.arange(300))
        
        if p2[0] > p2[7]:
            maxy = p2[0]
        else:
            maxy = p2[7]

        plt.xlim(xmin=0)
        plt.xlim(xmax=2100)
        #plt.ylim(ymax=maxy+50)

        plt.show()
        #pp.savefig()
        ax.clear()
        
        m = p2[3]
        b = p2[4]
        
        f.write("%9s \t %7f \t %6f \t%1d \n" %(pn,m,b,flag))
        
            # 0 = height of primary blue peak
            # 1 = location of primary blue peak
            # 2 = sigma of blue primary peak
            # 3 = linear term of dispersion
            # 4 = constant term of dispersion
            # 5 = scale factor for sigma
            # 6 = scale factor for Blue peaks
            # 7 = height of IR peak
        
        fp.write("%9s \t %1d \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \t %6f \n" %(pn,flag,parinfo[0]['limits'][0],p2[0],parinfo[0]['limits'][1],parinfo[1]['limits'][0],p2[1],parinfo[1]['limits'][1],parinfo[2]['limits'][0],p2[2],parinfo[2]['limits'][1],parinfo[3]['limits'][0],p2[3],parinfo[3]['limits'][1],parinfo[4]['limits'][0],p2[4],parinfo[4]['limits'][1],parinfo[5]['limits'][0],p2[5],parinfo[5]['limits'][1],parinfo[6]['limits'][0],p2[6],parinfo[6]['limits'][1],parinfo[7]['limits'][0],p2[7],parinfo[7]['limits'][1]))
        
        f.close()
        fp.close()

#pp.close()
#print "Skipped ", skipped
#print skiplist
#fskip = open("junk.cfg", "w")
#for pixel in skiplist:
    #fskip.write("pixel="+str(pixel)+"\n")
#f=open(outfile, 'a')
#f.write("dir="+str(dir)+"\n")
#f.write("20110804/obs_20110804-223423.h5")
#fskip.write("dir="+str(dir)+"\n")
#fskip.write("20110804/obs_20110804-223423.h5")
#fskip.close()
h5file.close()

