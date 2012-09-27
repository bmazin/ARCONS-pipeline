import numpy as np
from tables import *
import time
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
import scipy.signal as signal
from scipy import optimize
import scipy.stats as stats


# Define the various classes and functions needed for the beam mapping
# Define BeamMapper class - this allows you to change the location of a peak if there is a mistake
class BeamMapper():    
    # Initialize variables needed within the class
    def __init__(self,xtime,ytime):
        self.crx = np.zeros((2024,xtime))
        self.cry = np.zeros((2024,ytime))
        self.crx1 = np.zeros((2024,xtime))
        self.cry1 = np.zeros((2024,ytime))
        self.crx2 = np.zeros((2024,xtime))
        self.cry2 = np.zeros((2024,ytime))
        self.crx3 = np.zeros((2024,xtime))
        self.cry3 = np.zeros((2024,ytime))
        self.crx4 = np.zeros((2024,xtime))
        self.cry4 = np.zeros((2024,ytime))
        self.crx5 = np.zeros((2024,xtime))
        self.cry5 = np.zeros((2024,ytime))
        self.crx6 = np.zeros((2024,xtime))
        self.cry6 = np.zeros((2024,ytime))
        self.flag = np.zeros(2024)
        self.peakpos = np.zeros((2,2024))
    # Try to find a peak position by manually selecting an approximate peak location
    def on_click(self,event):
        # If x sweep plot (top plot) is clicked
        if(event.y > 250):
            self.xvals=np.arange(len(self.crx[pixelno][:]))
            self.xpeakguess=event.xdata
            self.xfitstart=max([self.xpeakguess-20,0])
            self.xfitend=min([self.xpeakguess+20,len(self.xvals)])
            params = fitgaussian(self.crx[pixelno][self.xfitstart:self.xfitend],self.xvals[self.xfitstart:self.xfitend])
            self.xfit = gaussian(params,self.xvals)
            self.peakpos[0,self.pixelno]=params[0]
        # If y sweep plot (bottom plot) is clicked
        else:
            self.yvals=np.arange(len(self.cry[pixelno][:]))
            self.ypeakguess=event.xdata
            self.yfitstart=max([self.ypeakguess-20,0])
            self.yfitend=min([self.ypeakguess+20,len(self.yvals)])
            params = fitgaussian(self.cry[pixelno][self.yfitstart:self.yfitend],self.yvals[self.yfitstart:self.yfitend])
            self.yfit = gaussian(params,self.yvals)
            self.peakpos[1,self.pixelno]=params[0]
    # Connect to plot
    def connect(self):
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)

# Define a standard Gaussian distribution function
def gaussian(pars, x):
    center, width, height, back = pars
    width = float(width)
    return back + height*np.exp(-(((center-x)/width)**2)/2)

# Define an error function between data and a Gaussian
def errorfunction(params, data, x):
    errorfunction = data - gaussian(params,x)
    return errorfunction

# Find an optimal Guassian fit for the data, return parameters of that Gaussian
def fitgaussian(data,x):
    params=(x.mean(),2.*(x[1]-x[0]),data.max(), 0.)
    p, success = optimize.leastsq(errorfunction, params, args=(data, x))
    return p

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


# Specify input/output directory and files
path = '/Users/paul/desktop/schoolwork/UCSB/Research/Sweep_2012-09-03/'
xsweep1 = 'Sweep_1/BM_X.h5'
ysweep1 = 'Sweep_1/BM_Y.h5'
xsweep2 = 'Sweep_2/BM_X.h5'
ysweep2 = 'Sweep_2/BM_Y.h5'
xsweep3 = 'Sweep_3/BM_X.h5'
ysweep3 = 'Sweep_3/BM_Y.h5'
xsweep4 = 'Sweep_4/BM_X.h5'
ysweep4 = 'Sweep_4/BM_Y.h5'
xsweep5 = 'Sweep_5/BM_X.h5'
ysweep5 = 'Sweep_5/BM_Y.h5'
xsweep6 = 'Sweep_6/BM_X.h5'
ysweep6 = 'Sweep_6/BM_Y.h5'

number_of_roaches = 8
roach_pixel_count = np.zeros(number_of_roaches)
for roachno in xrange(0,number_of_roaches):
    roach_pixel_count[roachno] = file_len(path + 'ps_freq%i.txt' %roachno)-1
print roach_pixel_count

# Load the input files
# X sweep data 1
h5file_x1 = openFile(path + xsweep1, mode = "r")
ts_x1 = int(h5file_x1.root.header.header.col('ut')[0])
exptime_x1 = int(h5file_x1.root.header.header.col('exptime')[0])
# Y sweep data 1
h5file_y1 = openFile(path + ysweep1, mode = "r")
ts_y1 = int(h5file_y1.root.header.header.col('ut')[0])
exptime_y1 = int(h5file_y1.root.header.header.col('exptime')[0])
# X sweep data 2
h5file_x2 = openFile(path + xsweep2, mode = "r")
ts_x2 = int(h5file_x2.root.header.header.col('ut')[0])
exptime_x2 = int(h5file_x2.root.header.header.col('exptime')[0])
# Y sweep data 2
h5file_y2 = openFile(path + ysweep2, mode = "r")
ts_y2 = int(h5file_y2.root.header.header.col('ut')[0])
exptime_y2 = int(h5file_y2.root.header.header.col('exptime')[0])
# X sweep data 3
h5file_x3 = openFile(path + xsweep3, mode = "r")
ts_x3 = int(h5file_x3.root.header.header.col('ut')[0])
exptime_x3 = int(h5file_x3.root.header.header.col('exptime')[0])
# Y sweep data 3
h5file_y3 = openFile(path + ysweep3, mode = "r")
ts_y3 = int(h5file_y3.root.header.header.col('ut')[0])
exptime_y3 = int(h5file_y3.root.header.header.col('exptime')[0])
# X sweep data 4
h5file_x4 = openFile(path + xsweep4, mode = "r")
#ts_x4 = int(h5file_x4.root.header.header.col('ut')[0])
ts_x4 = 1346722708
exptime_x4 = int(h5file_x4.root.header.header.col('exptime')[0])
# Y sweep data 4
h5file_y4 = openFile(path + ysweep4, mode = "r")
#ts_y4 = int(h5file_y4.root.header.header.col('ut')[0])
ts_y4 = 1346723385
exptime_y4 = int(h5file_y4.root.header.header.col('exptime')[0])
# X sweep data 5
h5file_x5 = openFile(path + xsweep5, mode = "r")
ts_x5 = int(h5file_x5.root.header.header.col('ut')[0])
exptime_x5 = int(h5file_x5.root.header.header.col('exptime')[0])
# Y sweep data 5
h5file_y5 = openFile(path + ysweep5, mode = "r")
ts_y5 = int(h5file_y5.root.header.header.col('ut')[0])
exptime_y5 = int(h5file_y5.root.header.header.col('exptime')[0])
# X sweep data 6
h5file_x6 = openFile(path + xsweep6, mode = "r")
ts_x6 = int(h5file_x6.root.header.header.col('ut')[0])
exptime_x6 = int(h5file_x6.root.header.header.col('exptime')[0])
# Y sweep data 6
h5file_y6 = openFile(path + ysweep6, mode = "r")
#ts_y6 = int(h5file_y6.root.header.header.col('ut')[0])
ts_y6 = 1346648400
exptime_y6 = int(h5file_y6.root.header.header.col('exptime')[0])


# Print start and sweep durations
print 'Start Time1 = ',ts_x1, ts_y1
print 'exptime_x =',exptime_x1,'and exptime_y =', exptime_y1
print '"A" = Accept, "D" = Delete, "Q" = Quit, "X" = X Only, "Y" = Y Only'

# Create a BeamMapper instance


# Go through each of the pixels (originally (0,2024))
for roachno in xrange(0,number_of_roaches):
    map = BeamMapper(exptime_x1,exptime_y1)
    f=open(path + 'r%i.pos' %roachno,'w')
    for pixelno in xrange(0,int(roach_pixel_count[roachno])):
        map.pixelno=pixelno
        map.flag[pixelno] = pixelno
    
        # Store the x data into crx  
        pn1 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_x1)
        pn2 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_x2)
        pn3 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_x3)
        pn4 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_x4)
        pn5 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_x5)
        pn6 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_x6)
        try:
            data1 = h5file_x1.root._f_getChild(pn1).read()
            data2 = h5file_x2.root._f_getChild(pn2).read()
            data3 = h5file_x3.root._f_getChild(pn3).read()
            data4 = h5file_x4.root._f_getChild(pn4).read()
            data5 = h5file_x5.root._f_getChild(pn5).read()
            data6 = h5file_x6.root._f_getChild(pn6).read()
            for j in xrange(0,exptime_x1):
                map.crx[pixelno][j] = np.median([len(data1[j]), len(data2[j]), len(data3[j]), len(data4[j]), len(data5[j]), len(data6[j+12])])
                map.crx1[pixelno][j] = len(data1[j])
                map.crx2[pixelno][j] = len(data2[j])
                map.crx3[pixelno][j] = len(data3[j])
                map.crx4[pixelno][j] = len(data4[j])
                map.crx5[pixelno][j] = len(data5[j])
                map.crx6[pixelno][j] = len(data6[j+12])
                
        except:
            pass
        
        # Store the y data into cry
        pn1 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_y1)
        pn2 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_y2)
        pn3 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_y3)
        pn4 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_y4)
        pn5 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_y5)
        pn6 = '/r%d/p%d/t%d' % ( roachno ,pixelno, ts_y6)
        try:
            data1 = h5file_y1.root._f_getChild(pn1).read()
            data2 = h5file_y2.root._f_getChild(pn2).read()
            data3 = h5file_y3.root._f_getChild(pn3).read()
            data4 = h5file_y4.root._f_getChild(pn4).read()
            data5 = h5file_y5.root._f_getChild(pn5).read()
            data6 = h5file_y6.root._f_getChild(pn6).read()
            for j in xrange(exptime_y1):
                map.cry[pixelno][j] = np.median([len(data1[j]), len(data2[j]), len(data3[j]), len(data4[j]), len(data5[j]), len(data6[j-29])])
                map.cry1[pixelno][j] = len(data1[j])
                map.cry2[pixelno][j] = len(data2[j])
                map.cry3[pixelno][j] = len(data3[j])
                map.cry4[pixelno][j] = len(data4[j])
                map.cry5[pixelno][j] = len(data5[j])
                map.cry6[pixelno][j] = len(data6[j-29])
                
        except:
            pass

        print 'roach %d, pixel %d' % (roachno, pixelno)
        map.fig = plt.figure()

        # do fit of x-position
        map.ax = map.fig.add_subplot(211)
        map.xvals=np.arange(len(map.crx[pixelno][:]))
        plt.title(str(pixelno))
        plt.plot(map.xvals,map.crx[pixelno][:])       
        map.xpeakguess=np.where(map.crx[pixelno][:] == map.crx[pixelno][:].max())[0][0]
        map.xfitstart=max([map.xpeakguess-20,0])
        map.xfitend=min([map.xpeakguess+20,len(map.xvals)])
        params_x = fitgaussian(map.crx[pixelno][map.xfitstart:map.xfitend],map.xvals[map.xfitstart:map.xfitend])
        print 'x: [center, width, height, back] =', params_x
        map.xfit = gaussian(params_x,map.xvals)
        plt.plot(map.xvals, map.xfit)
        plt.plot(map.xvals, map.crx1[pixelno][:],alpha = .2)
        plt.plot(map.xvals, map.crx2[pixelno][:],alpha = .2)
        plt.plot(map.xvals, map.crx3[pixelno][:],alpha = .2)
        plt.plot(map.xvals, map.crx4[pixelno][:],alpha = .2)
        plt.plot(map.xvals, map.crx5[pixelno][:],alpha = .2)
        plt.plot(map.xvals, map.crx6[pixelno][:],alpha = .2)
        map.ax.set_ylim(params_x[3]-30,params_x[2]+75)
        map.peakpos[0,pixelno]=params_x[0]
            
        # do fit of y-position
        map.ax = map.fig.add_subplot(212)
        map.yvals=np.arange(len(map.cry[pixelno][:]))
        plt.title(str(pixelno+1))
        plt.plot(map.yvals,map.cry[pixelno][:])     
        map.ypeakguess=np.where(map.cry[pixelno][:] == map.cry[pixelno][:].max())[0][0]
        map.yfitstart=max([map.ypeakguess-20,0])
        map.yfitend=min([map.ypeakguess+20,len(map.yvals)])
        params_y = fitgaussian(map.cry[pixelno][map.yfitstart:map.yfitend],map.yvals[map.yfitstart:map.yfitend])
        print 'y: [center, width, height, back] =', params_y
        map.yfit = gaussian(params_y,map.yvals)
        plt.plot(map.yvals, map.yfit)
        plt.plot(map.yvals,map.cry1[pixelno][:],alpha = .2)
        plt.plot(map.yvals,map.cry2[pixelno][:],alpha = .2)
        plt.plot(map.yvals,map.cry3[pixelno][:],alpha = .2)
        plt.plot(map.yvals,map.cry4[pixelno][:],alpha = .2)
        plt.plot(map.yvals,map.cry5[pixelno][:],alpha = .2)
        plt.plot(map.yvals,map.cry6[pixelno][:],alpha = .2)
        map.ax.set_ylim(params_y[3]-30,params_y[2]+75)
        map.peakpos[1,pixelno]=params_y[0]

        # AUTOMATICALLY DELETE OBVIOUSLY BAD PIXELS
        
#        if((params_x[0] == 9.5 and params_y[0] == 9.5) or (params_x[1] >= 100 and params_y[1] >= 100)):
#            map.peakpos[0,pixelno]=0
#            map.peakpos[1,pixelno]=0
#            map.flag[pixelno] = -1*pixelno
#            plt.close()
#            f=open(path+'Sweep_1/r%i.pos' %roachno,'a')
#            f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+'\n')
#            f.close()
#        else:
        map.connect()
        while True:
            # Accept pixel with 'a'
            response=raw_input('>')
            if(response == 'a'):
                print 'Response: ', response
                plt.close()
                f=open(path+'r%i.pos' %roachno,'a')
                f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+'\n')
                f.close()
                break
            # Delete pixel with 'd'
            elif(response == 'd'):
                print 'Response: ', response
                map.peakpos[0,pixelno]=0
                map.peakpos[1,pixelno]=0
                map.flag[pixelno] = -1*pixelno
                plt.close()
                f=open(path+'r%i.pos' %roachno,'a')
                f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+'\n')
                f.close()
                break
            # Control whether only the x or y pixel is selected
            # x is okay
            elif(response == 'x'):
                print 'Response: ', response
                map.peakpos[1,pixelno]=0
                map.flag[pixelno] = -1*pixelno
                plt.close()
                f=open(path+'r%i.pos' %roachno,'a')
                f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+'\n')
                f.close()
                break
            # y is okay
            elif(response == 'y'):
                print 'Response: ', response
                map.peakpos[0,pixelno]=0
                map.flag[pixelno] = -1*pixelno
                plt.close()
                f=open(path+'r%i.pos' %roachno,'a')
                f.write(str(map.peakpos[0,pixelno])+'\t'+str(map.peakpos[1,pixelno])+'\t'+str(int(map.flag[pixelno]))+'\n')
                f.close()
                break
            # Quit program with 'q'
            elif(response == 'q'):
                sys.exit(0)
            plt.show()
            
h5file_x.close()
h5file_y.close()
f.close()
