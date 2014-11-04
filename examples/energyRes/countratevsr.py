'''
Author: Danica Marsden                          Date: June 4, 2013

'''

import sys, os
import time
import tables
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import mpfit
import smooth
from utils import *
from tables import *
from fitFunctions import *
from matplotlib.backends.backend_pdf import PdfPages

## Constants:

h = 4.1357e-15   # eV s    
c = 3.e8        # m/s
ang2m = 1.e-10
fwhm2sig = 2.35

n_rows = 46
n_cols = 44
n_roaches = 8
non_alloc_pix = 250

n_bits_base = 20
n_bits_peak = 32
n_bits_parab = 44 
n_bits_long = 12
n_bits_time = 20

#cr_high = 2499.
cr_high = 3000.
dead_time = 1.e-4

big_step = 10       
bin_smooth = 60.


pulseMask = int(n_bits_long*'1',2)   # bitmask of 12 ones
timeMask = int(n_bits_time*'1',2)

#time_bins_start = [110.,235.,360.,490.,615.,745.,870.,995.,1120.,1250.,1375.,1500.,1630.,1755.,1880.,2005.,2135.,2260.,2390.,2515.,2640.,2765.,2895.,3020.,3145.,3270.]
#time_bin_width = 45.
#time_bins_start = [108.,195.,281.,367.,453.,540.,626.,712.,800.,886.,973.,1058.,1145.,1231.,1317.,1403.,1490.,1574.]
#time_bin_width = 18.
time_bins_start = [74.,201.,327.,454.,581.,707.,833.,960.,1087.,1214.,1340.,1467.,1594.,1721.,1846.]
time_bin_width = 56.


print 'Start Time: ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())


## Read in the scaled count rates we should be seeing:

incounts = []

#powermeter = open('correctedCountRateStair550nm.dat', 'r')
#powermeter = open('correctedCountRate550nm.dat', 'r')
powermeter = open('obs_20130612-003423_CorrectedCount.txt', 'r')
        
for line in powermeter:
    if '#' in line:
        continue

    incounts.append(float(line.strip()))
    
incounts = np.array(incounts)
powermeter.close()


## For each pixel, get the corresponding count rates, and R,
## and make the plot cts/s vs cts/s and R


#h5file = '/Scratch/linearityTestData/obs_20130531-002421.h5'
#h5file = '/Scratch/linearityTestData/obs_20130530-223349.h5'
obsFileName = '/Scratch/linearityTestData/obs_20130612-003423.h5'



obs = ObsFile(obsFileName)
# open the file
data = tables.openFile(h5file, mode='r')


beammap = data.root.beammap.beamimage.read()  # array of strings that translates
          # physical location [row_i,col_j] to the roach/pixel allocation in the
          # readout
    
#pp = PdfPages("countRatevsR_2.pdf")
#mpl.rcParams['font.size']=4
#n_plots_per_page = 20
#plotcounter = 0

#pixtouse = [[0,26],[0,30],[0,34],[1,28],[1,39],[2,25],[2,27],[2,31],[3,17],[3,26],[3,36],[4,28],[4,32],[5,32],[5,33],[6,22],[6,30],[7,24],[7,27],[7,28],[7,32],[12,26],[15,24],[16,23],[18,24],[19,22],[19,23],[19,24],[20,24],[23,24]]
pixtouse = [[0,26],[1,39],[2,25],[2,27],[2,31],[2,34],[2,39],[3,17],[3,26],[3,32],
            [3,36],[4,26],[4,32],[5,32],[6,22],[6,29],[6,30],[7,25],[7,27],[7,28],
            [7,32],[7,37],[8,30],[10,32],[10,42],[11,30],[12,29],[12,31],[13,29],[13,31],
            [14,27],[14,32],[15,25],[15,26],[15,30],[15,31],[15,32],[16,23],[16,24],[16,30],
            [17,23],[17,25],[17,26],[17,39],[18,24],[18,28],[19,23],[19,26],[19,28],[19,30]]  

avg_ct = np.zeros(len(incounts))
avg_ct_corr = np.zeros(len(incounts))
avg_r = np.zeros(len(incounts))
avg_r_err = np.zeros(len(incounts))
avg_counter = np.zeros(len(incounts))
for i in range(n_rows):
    for j in range(n_cols):

        pstring = beammap[i][j]

        if [i,j] in pixtouse:
            print i, j
        else:
            continue

        #print '\n----------------------------------------------\n', pstring, i, j
        #print i,j
        
        pixelNode = pstring.split('t')[0]

        roach = int(pixelNode.split('/')[1][1:])
        pixel = int(pixelNode.split('/')[2][1:])

        # cut non-allocated pixels
        if (roach==0) & (pixel==non_alloc_pix):
            print "Non allocated pixel"
            continue

        # lab testing - just roaches 0-3
        if (roach > 3):
            print "Non allocated pixel (roach not used)"
            continue

        # the non-sorted list of all photon packets belonging to pixel i,j
        photondata = np.concatenate(data.root._f_getChild(pstring)[:])

        # flag dead pixels
        if len(photondata)==0:
            print "Dead pixel"
            continue

        # for each second, get the packets
        timeinsec = []
        count_rate = []
        first_one = 0
        for sec in range(len(data.getNode(pstring))):
            
            packets = data.getNode(pstring)[sec]
            
            microtimes = packets & timeMask
            peak_phase_sub = packets >> n_bits_peak & pulseMask
            parab_phase_sub = packets >> n_bits_parab & pulseMask
            baseline_sub = packets >> n_bits_base & pulseMask

            times = float(sec) + microtimes*(1.e-6)

            # count rate
            cr = len(packets)
            count_rate.append(float(cr))    # per second
            timeinsec.append(float(sec))

            # put it all together
            if first_one == 0:
                peak_phase = peak_phase_sub
                parab_phase = parab_phase_sub
                baseline = baseline_sub
                timestamp = times
                first_one = 1
            else:
                peak_phase = np.concatenate((peak_phase,peak_phase_sub))
                parab_phase = np.concatenate((parab_phase,parab_phase_sub))
                baseline = np.concatenate((baseline,baseline_sub))
                timestamp = np.concatenate((timestamp,times))

        timestamp = np.array(timestamp,dtype=float)
        peak_phase = np.array(peak_phase, dtype=float)
        parab_phase = np.array(parab_phase, dtype=float)
        baseline = np.array(baseline, dtype=float)
        photon_ind = np.arange(baseline.size, dtype=float)
        
        timeinsec = np.array(timeinsec, dtype=float)
        count_rate = np.array(count_rate, dtype=float)
        

        # subtract out the baseline
        parab_phase_orig = parab_phase
        #parab_phase -= baseline            
        peak_subt = peak_phase - baseline

        ## data sliced out where pixels go hot
        #fig = plt.figure()
        #plt.plot(timeinsec, count_rate)
        #plt.ylabel('Count Rate [#/s]')
        #plt.xlabel('Time [s]')
        #plt.show()
        #plt.clf()


        #fig = plt.figure()
        
        #ax1 = fig.add_subplot(2,1,1)
        #plt.plot(timestamp, baseline, '.', markersize=2)
        #plt.ylabel('Baseline Amplitude')
        #plt.xlabel('Time [s]')
        #plt.xlim(0.,800.)
        #plt.ylim(800.,2100.)

        #ax2 = fig.add_subplot(3,1,2)
        #plt.plot(timestamp, peak_phase, '.', markersize=1)
        #plt.ylabel('Peak')
        #plt.xlabel('Time [s]')
        #plt.xlim(0.,800.)
        #plt.ylim(800.,2050.)        
        
        #ax3 = fig.add_subplot(2,2,3)
        #plt.plot(timestamp, parab_phase_orig, '.', markersize=1)
        #plt.ylabel('Parab')
        #plt.xlabel('Time [s]')
        #plt.xlim(0.,800.)
        #plt.ylim(800.,2050.)

        #ax4 = fig.add_subplot(2,1,2)
        #plt.plot(timestamp, peak_subt, '.', markersize=2)
        #plt.ylabel('Peak-baseline')
        #plt.xlabel('Time [s]')
        #plt.xlim(0.,800.)
        #plt.ylim(-1100.,0.)
        
        #plt.show()
        #plt.clf()
        
        timestamp_cut = []
        baseline_cut = []
        peak_phase_cut = []
        peak_subt_cut = []
        for qq in range(len(timestamp)-2):
            if (timestamp[qq+1] - timestamp[qq]) > 0.001:
                timestamp_cut.append(timestamp[qq+1])
                baseline_cut.append(baseline[qq+1])
                peak_phase_cut.append(peak_phase[qq+1])
                peak_subt_cut.append(peak_subt[qq+1])

        #fig = plt.figure()
        
        #ax1 = fig.add_subplot(3,1,1)
        #plt.plot(timestamp_cut, baseline_cut, '.', markersize=1)
        #plt.ylabel('Baseline Amplitude')
        #plt.xlabel('Time [s]')
        #plt.title('After cut dt=200us')
        #plt.xlim(0.,600.)
        #plt.ylim(1500.,2200.)

        #ax2 = fig.add_subplot(3,1,2)
        #plt.plot(timestamp_cut, peak_phase_cut, '.', markersize=1)
        #plt.ylabel('Peak')
        #plt.xlabel('Time [s]')
        #plt.xlim(0.,600.)
        #plt.ylim(800.,2050.)        
        
        #ax3 = fig.add_subplot(3,1,3)
        #plt.plot(timestamp_cut, peak_subt_cut, '.', markersize=1)
        #plt.ylabel('Peak-baseline')
        #plt.xlabel('Time [s]')
        #plt.xlim(0.,600.)
        #plt.ylim(-1100.,0.)
        
        #plt.show()
        #plt.clf()


        timestamp = timestamp_cut
        timestamp = np.array(timestamp,dtype=float)
        baseline = baseline_cut
        baseline = np.array(baseline,dtype=float)
        peak_phase = peak_phase_cut
        peak_phase = np.array(peak_phase,dtype=float)
        peak_subt = peak_subt_cut
        peak_subt = np.array(peak_subt,dtype=float)

        pixel_counts = []
        resolution = []
        resolution_error = []
        ticker = 0
        which_ones = []
        labels = []
        plotcounter = 0

        baseline_avgs = []
        dtimes = []
        dbase = []
        dsignal = []
        ## come up with an average count rate in each interval, cutting hot bits
        for timebin in time_bins_start:

            #print ticker
            ticker += 1

            condition = (timeinsec >= timebin)
            
            counts_in_bin = count_rate[condition][0:time_bin_width]
            times_in_bin = timeinsec[condition][0:time_bin_width]

            ## cut hot bits
            counts_in_bin = counts_in_bin[np.where(counts_in_bin < cr_high)[0]]
            hot_times = times_in_bin[np.where(counts_in_bin > cr_high)[0]]

            ## average count rate
            avg_countrate = np.sum(counts_in_bin)/len(counts_in_bin)
            pixel_counts.append(avg_countrate)
            if avg_countrate > 500.:
                low_counts_val = 20.
            else:
                low_counts_val = 2.

            condition2 = np.where((timestamp > timebin) & (timestamp < timebin+time_bin_width))[0]
            
            timestamps_bin = timestamp[condition2]
            phases_in_bin = parab_phase[condition2]
            peak_in_bin = peak_phase[condition2]
            parab_phase_orig_bin = parab_phase_orig[condition2]
            baseline_bin = baseline[condition2]


            #print 'Starting dt loop', ticker, ' at ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
            baseline_avgs.append(np.mean(np.array(baseline_bin)))
            dtimes.append(timestamps_bin[1] - timestamps_bin[0])
            dbase.append(abs(baseline_bin[0]-baseline_avgs[-1]))
            for q in range(len(timestamps_bin)-2):
                #dtimes.append(min((timestamps_bin[q+2]-timestamps_bin[q+1]), (timestamps_bin[q+1]-timestamps_bin[q])))
                dtimes.append(timestamps_bin[q+1] - timestamps_bin[q])
                dbase.append(abs(baseline_bin[q+1]-baseline_avgs[-1]))
            #dtimes.append(timestamps_bin[len(timestamps_bin)-1] - timestamps_bin[len(timestamps_bin)-2])
            #dbase.append(abs(baseline_bin[len(timestamps_bin)-1]-baseline_avgs[-1]))       
            #print 'Ended dt loop', ticker, ' at ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())


            if hot_times != []:
                print hot_times
                good_times = []
                good_phases = []
                for t in range(len(timestamps_bin)):
                    if np.floor(timestamps_bin[t]) in hot_times:
                        continue
                    good_times.append(timestamps_bin[t])
                    good_phases.append(phases_in_bin[t])
                good_times = np.array(good_times)
                good_phases = np.array(good_phases)
            else:
                good_times = timestamps_bin
                good_phases = phases_in_bin



            ## For baseline testing:
            #good_phases = phases_in_bin - baseline_bin
            #good_phases = peak_in_bin
            good_phases = peak_in_bin - baseline_bin
            #good_phases = baseline_bin
            good_times = timestamps_bin

        

            ## histogram phases in each time interval

            try:
                rangex = max(good_phases)-min(good_phases)
                nbins = int(np.ceil(np.abs(rangex)))
                n_inbin, phasebins = np.histogram(good_phases, nbins)
            except:
                print 'histogram not working'
                pixel_counts.pop()
                continue

            binwidth = phasebins[1]-phasebins[0]                    # 1.0
            phasebins += binwidth/2.0
            phasebins = phasebins[0:len(phasebins)-1]
            n_inbin = np.array(n_inbin)                            
            phasebins = np.array(phasebins)

            ## smooth

            window = 'hanning'
            windowsize = int(bin_smooth)    ## needs to be even to work -- fix in smooth?
            try:
                parab_smooth = smooth.smooth(n_inbin, windowsize, window)
                smoothed_data = np.array(parab_smooth, dtype=float)
                start_indsm = (parab_smooth > low_counts_val).nonzero()[0][0]
                end_indsm = (parab_smooth > low_counts_val).nonzero()[0][-1]
            except:
                print 'smooth didnt work'
                pixel_counts.pop()
                continue 

            ## find extreema:

            coarse_data_len = np.floor(len(smoothed_data)/big_step)
            coarse_data = np.zeros(coarse_data_len)
            coarse_x = np.zeros(coarse_data_len)
       
            for k in range(len(coarse_data)):
                coarse_data[k] = smoothed_data[k*big_step]
                coarse_x[k] = np.array(phasebins)[k*big_step]

            try:
                start_ind = (coarse_data > low_counts_val).nonzero()[0][0]
                end_ind = (coarse_data > low_counts_val).nonzero()[0][-1]
            except:
                print 'index ends didnt work'
                pixel_counts.pop()
                continue

                
            coarse_data = coarse_data[start_ind:end_ind]
            coarse_x = coarse_x[start_ind:end_ind]
            
            gradients = np.diff(coarse_data)

            maxima_num = 0
            minima_num = 0
            max_locations = []
            min_locations = []
            max_vals = []
            min_vals = []
            count = 0
            for grad in gradients[:-1]:
                count += 1

                if ((grad > 0) & (gradients[count] < 0) & (grad != gradients[count]) & (coarse_data[count] > 1)):
                    maxima_num += 1
                    max_locations.append(coarse_x[count])
                    max_vals.append(coarse_data[count])

                if ((grad < 0) & (gradients[count] > 0) & (grad != gradients[count])):
                    minima_num += 1
                    min_locations.append(coarse_x[count])
                    min_vals.append(coarse_data[count])


            #print 'maxima: ', max_locations, 'max vals: ', max_vals
            #print 'minima: ', min_locations, 'min vals: ', min_vals

            min_locations = np.array(min_locations)
            max_locations = np.array(max_locations)

            try:
                ind_right = np.max(np.where(phasebins < min_locations[-1]))
                ind_left = np.min(np.where(phasebins > coarse_x[0]))

                #ind_left = np.min(np.where(n_inbin > (max_vals[0]/2.)))
                #ind_right = np.max(np.where(n_inbin > (max_vals[0]/2.)))
            except:
                print 'now this one didnt work'
                pixel_counts.pop()
                continue
            
            ## fit a Gaussian

            amplitude = max(n_inbin[ind_left:ind_right])
            #x_offset = np.median(phasebins[ind_left:ind_right])
            #sigma = 50.
            x_offset = max_locations[0]
            sigma = np.abs(phasebins[ind_left] - phasebins[ind_right])/2.
            #print 'sigma: ', sigma
            sigma_orig = sigma
            y_offset = 1.e-8
        
            params=[sigma, x_offset, amplitude, y_offset]  # First guess at fit params
            errs = np.sqrt(n_inbin[ind_left:ind_right])                         # Poisson counts 
            errs[np.where(errs == 0.)] = 1.
            quiet = True

            parinfo = [ {'n':0,'value':params[0],'limits':[1., 500.], 'limited':[True,True],'fixed':False,'parname':"Sigma",'error':0},
               {'n':1,'value':params[1],'limits':[x_offset-sigma*2, x_offset+sigma*2],'limited':[True,True],'fixed':False,'parname':"x offset",'error':0},
               {'n':2,'value':params[2],'limits':[1., 2.*amplitude],'limited':[True,True],'fixed':False,'parname':"Amplitude",'error':0},
               {'n':3,'value':params[3],'limited':[False,False],'fixed':False,'parname':"y_offset",'error':0}]

            fa = {'x':phasebins[ind_left:ind_right],'y':n_inbin[ind_left:ind_right],'err':errs}

            m = mpfit.mpfit(gaussian, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
            if (quiet == False):
                print m.status, m.errmsg
            
            mpp = m.params                                #The fit params
            mpperr = m.perror
            
            for k,p in enumerate(mpp):
                parinfo[k]['value'] = p
                #print parinfo[k]['parname'],p," +/- ",mpperr[j]
                if k==0: sigma = p
                if k==1: x_offset = p
                if k==2: amplitude = p
                if k==3: y_offset = p

            gaussfit = y_offset + amplitude * np.exp( - (pow(( phasebins[ind_left:ind_right] - x_offset),2) / ( 2. * pow(sigma,2))))
            #print sigma

            if plotcounter==0:
                fig3=plt.figure()
                ax2=fig3.add_subplot(111)
                colormap = mpl.cm.gist_ncar
                plt.gca().set_color_cycle([colormap(z) for z in np.linspace(0, 0.9, 15)])
                ax2.plot(phasebins, n_inbin/float(np.max(n_inbin)))
                labels.append('count rate: %i, sigma: %1.f' %(avg_countrate, sigma))
                plt.title('[%i,%i]' %(i,j))
                plotcounter+=1
            else:
                ax2.plot(phasebins, n_inbin/float(np.max(n_inbin)))#,color = plt.get_cmap('jet')(float(iterator)/(n_roaches)))
                labels.append('count rate: %i, sigma: %1.f' %(avg_countrate, sigma))

            #print avg_countrate, sigma, x_offset, x_offset/sigma
            #fig = plt.figure()
            #plt.plot(phasebins, n_inbin)
            #plt.plot(phasebins, parab_smooth, 'y')
            #plt.plot(max_locations, max_vals, 'ko', markersize=10)
            #plt.plot(min_locations, min_vals, 'mo', markersize=10)
            #plt.plot(phasebins[ind_left:ind_right], gaussfit,'r')
            #plt.show()
            #plt.clf()
            
            ## get R

            resolution.append(np.abs(x_offset)/(2.355*sigma))
            resolution_error.append(np.abs(x_offset)/(2.355*sigma) * np.sqrt((mpperr[0]/sigma)**2 + (mpperr[1]/x_offset)**2))
            #print np.abs(x_offset)/(2.355*sigma), incounts[ticker-1], avg_countrate

            which_ones.append(ticker-1)

        plt.legend(labels, loc='upper left')
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()  
        llines = leg.get_lines() 
        plt.setp(ltext, fontsize='small')    # the legend text fontsize
        plt.setp(llines, linewidth=1.5)
        #plt.xlim(-1000.,0.)
        #plt.xlabel('Phase-Base')
        #plt.savefig(figdir+outfile+"_R_Estimates.png")
        plt.show()
        plt.clf()

        print len(time_bins_start), len(baseline_avgs), len(dtimes), len(dbase), len(timestamp)
        print timestamp[0:10]
        print dtimes[0:10]
        #dtimes_sorted = dtimes
        #dtimes_sorted.sort()
        #print dtimes_sorted[0:30]

        print baseline_avgs

        #fig = plt.figure()
        #plt.plot(dtimes, dbase, '.', markersize=2)
        #plt.xlim(-0.001,0.08)
        #plt.title('[%i,%i]' %(i,j))
        #plt.xlabel('Time between pulse and one before [s]')
        #plt.ylabel('Baseline - Average')
        #plt.show()
        #plt.clf()

        
        pixel_counts = np.array(pixel_counts)
        corr_count_rate = pixel_counts + pixel_counts*(dead_time*pixel_counts) + pixel_counts*(dead_time*pixel_counts)**2

        resolution = np.array(resolution)
        resolution_error = np.array(resolution_error)

        length = len(pixel_counts)
        if (length <5):
            continue
        
        incounts_sub = incounts[which_ones[:]]

        avg_ct[which_ones[:]] += pixel_counts
        avg_ct_corr[which_ones[:]] += corr_count_rate
        avg_r[which_ones[:]] += resolution * 1.375   ## 550nm to 400nm
        avg_r_err[which_ones[:]] += resolution_error
        avg_counter[which_ones[:]] += np.ones((len(which_ones)))

        ## order them
        ord = incounts_sub.argsort()
        incounts_sub = incounts_sub[ord]
        pixel_counts = pixel_counts[ord]
        corr_count_rate = corr_count_rate[ord]
        resolution = resolution[ord]


        ## Plot counts/s in vs. counts/s at detector and R


        #if (plotcounter % n_plots_per_page == 0):
        #    fig6 = plt.figure(figsize=(8.25, 10), dpi=100)
            
        ##plt.subplot(5,4,plotcounter%20+1)
        ##ax1 = fig6.add_subplot(111)
        #ax1 = fig6.add_subplot(5,4,plotcounter%20+1)

        #fig6 = plt.figure()
        #ax1 = fig6.add_subplot(111)
        
        #ax1.plot(incounts_sub,pixel_counts)
        #ax1.plot(incounts_sub,corr_count_rate,'g--')

        #ax2 = ax1.twinx()
        #ax2.plot(incounts_sub,resolution, 'r')
            
        #titlestring = str(i)+', '+str(j)+' '+pstring
        #plt.title(titlestring)
        #plt.show()
        #plt.clf()
        
        #if (((plotcounter +1) % n_plots_per_page == 0) or ((i+1 == n_rows) and (j+1 == n_cols))):
        #    pp.savefig(fig6)
            
        #plotcounter+=1
            
       
        ## make final pdf plot over all pixels that are "good"
        ## also counts/s out and R as a function of dist from center
        ## email out

#pp.close()
#print plotcounter, ' good pixels'

avg_ct /= avg_counter
avg_ct_corr /= avg_counter
avg_r /= avg_counter
avg_r_err /= avg_counter 

#n_avg = 14
n_avg = 15

for k in range(n_avg):
    print incounts[k], avg_ct[k], avg_ct_corr[k], avg_r[k], avg_counter[k], avg_r_err[k]


fig7 = plt.figure()
ax1 = fig7.add_subplot(111)
        
ax1.plot(incounts[0:n_avg],avg_ct[0:n_avg])
ax1.plot(incounts[0:n_avg],avg_ct_corr[0:n_avg],'g--')
ax1.set_ylabel('Measured counts [photons/s]')

ax2 = ax1.twinx()
ax2.plot(incounts[0:n_avg],avg_r[0:n_avg], 'r')
ax2.set_ylabel('Resolution at 400nm')
            
ax1.set_xlabel('Input counts [photons/s]')
plt.title('Average over 30 pixels')
plt.show()
plt.clf()

print 'Finish Time: ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

data.close()

