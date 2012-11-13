'''
Author: Danica Marsden                         Date: August 15, 2012

Reads in a wavelength calibration file (blue, red and IR laser light
) that has a beam map solution, and makes a histogram of photon phase
amplitudes for each pixel. Bad pixels are flagged.  Gaussians are fit
to the histograms.

The peaks are fit with a polynomial and that gives a phase amplitude
<-> wavelength correspondance.

Returns the polynomial fit parameters for each pixel, the error in
Angstroms (the width of the Gaussian), and the start and stop wave-
lengths, as well as the array of flagged pixels.

To run inside another script:

paramFile = sys.argv[1]
wavelengthCal(paramFile)

where paramFile has the same parameters as e.g. waveCal.dict

Depends on readDict.py, mpfit.py, smooth.py, fitFunctions.py

see README.txt for more info
'''
#!/bin/env python

import sys, os
import time
import tables
import numpy as np
import matplotlib as mpl
import matplotlib.pylab as plt
import mpfit
import smooth
import readDict
from utils import *
from tables import *
from fitFunctions import *
from matplotlib.backends.backend_pdf import PdfPages


class WaveCalSoln(IsDescription):
    roach = tables.UInt16Col()                 # ROACH board number
    pixelnum = tables.UInt16Col()              # pixel number on the roach
    pixelrow = tables.UInt16Col()              # physical x location - from beam map
    pixelcol = tables.UInt16Col()              # physical y location 
    polyfit = tables.Float64Col(3)             # polynomial to convert from phase amplitude to wavelength,
                                               #    double precision
    sigma = tables.Float64Col()                # 1 sigma (Gaussian width) in eV, for blue peak
    solnrange = tables.Float32Col(2)           # start and stop wavelengths for the fit in Angstroms
    wave_flag = tables.UInt16Col()             # flag to indicate if pixel is good (0), dead (1) or failed
                                               #    during wave cal fitting (2)

class DriftObj(IsDescription):
    pixelrow = tables.UInt16Col()              # physical x location - from beam map
    pixelcol = tables.UInt16Col()              # physical y location 
    gaussparams = tables.Float64Col(6)         # polynomial to convert from phase amplitude to wavelength,
                                               #    double precision
    
def failure(row, rarray, flagnum):
    row['wave_flag'] = flagnum
    row['polyfit'] = np.array([-1.,-1., -1.])
    row['sigma'] = -1.
    row['solnrange'] = np.array([-1.,-1.])
    rarray.append(0.)

    
def wavelengthCal(paramFile): 


    params = readDict.readDict(paramFile)
    params.readFromFile(paramFile)

    infile = params['calfileslist']
    files2run = params['files2run']
    outdir = params['outdir']

    cal_files = []

    file_list = open(outdir+infile)

    if 'cal_' in files2run:
        cal_files.append(files2run)
    else:
        if '201209' in files2run:
            while 1:
                line = file_list.readline()
                if files2run in line:
                    cal_files.append(line)
        else:
            while 1:
                line = file_list.readline()
                cal_files.append(line)
    

    for k in range(len(cal_files)):

        print 'Starting file ', cal_files[k], ' at time: ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

        datedir = cal_files[k].split('/')[3]+'/'
        try:
            os.mkdir(outdir+datedir)
        except:
            pass
    
        data = tables.openFile(cal_files[k],mode='r')

        n_rows = params['n_rows']
        n_cols = params['n_cols']
    
        xarray = []
        yarray = []
        rarray = []

        # open outfile to write parameters to
        outfile = 'calsol' + cal_files[k].split('cal')[1]
        try:
            h5out = tables.openFile(outdir+datedir+outfile, mode='w')
        except:
            pass
        root = h5out.root

        # create a group that branches from the root node, save cal data table in this group
        calgroup = h5out.createGroup(root, 'wavecal', 'Table of calibration parameters for each pixel')

        # create a table instance under group calgroup with node name "calsoln"
        # the WaveCalSoln class declared before is the description parameter to define the columns of the table
        caltable = h5out.createTable(calgroup, 'calsoln', WaveCalSoln, title='Wavelength Cal Table')

        # pointer to the table row instance:
        row = caltable.row
    
        # Do the same for the file to monitor pixel drifts
        try:
            os.mkdir(outdir+datedir+params['driftdir'])
        except:
            pass
        try:
            h5outdrift = tables.openFile(outdir+datedir+params['driftdir']+outfile.split('.')[0]+'_drift.h5', mode='w')
        except:
            pass
        driftroot = h5outdrift.root
        driftgroup = h5outdrift.createGroup(driftroot, 'params_drift', 'Table of parameters for drift study')
        drifttable = h5outdrift.createTable(driftgroup, 'driftparams', DriftObj, title='Drift Params Table')
        driftrow = drifttable.row


        # get the beam map, an array of strings that translates physical location [row_i,col_j] to the roach/pixel in the readout
        beammap = data.root.beammap.beamimage.read()  


        # plot of fits for each "good" pixel
        try:
            os.mkdir(outdir+datedir+params['figdir'])
        except:
            pass
        
        pp = PdfPages(outdir+datedir+params['figdir']+outfile.split('.')[0]+"_fits.pdf")
        mpl.rcParams['font.size'] = 4
        n_plots_per_page = 20
        plotcounter = 0


        # loop through each pixel and get its solution, populate output WaveCalSoln and append the row.
        for i in range(n_rows):
            for j in range(n_cols):

                pstring = beammap[i][j]
                pixelNode = pstring.split('t')[0]
    
                xarray.append(i)
                yarray.append(j)
                roach = int(pixelNode.split('/')[1][1:])
                pixel = int(pixelNode.split('/')[2][1:])

                row['pixelrow'] = i
                row['pixelcol'] = j
                row['roach'] = roach
                row['pixelnum'] = pixel
                row['wave_flag'] = 0                      # 0 until flagged (1 for dead, 2 for bad data/fit)

                # the non-sorted list of all photon packets belonging to pixel i,j
                photondata = np.concatenate(data.root._f_getChild(pstring)[:])

                # flag dead pixels
                if len(photondata)==0:
                    failure(row, rarray, 1)
                    continue

                # for each second, get the packets
                timeinsec = []
                count_rate = []
                first_one = 0
                n_cut = 0
                for sec in range(len(data.getNode(pstring))):
            
                    packets = data.getNode(pstring)[sec]
            
                    microtimes = packets & params['timeMask']
                    parab_phase_sub = packets >> params['n_bits_parab'] & params['pulseMask']
                    baseline_sub = packets >> params['n_bits_base'] & params['pulseMask']

                    times = float(sec) + microtimes*(1.e-6)

                    # cut on count rate - data sliced out where pixels go hot
                    cr = len(packets)
                    if ((cr > params['cr_low'])*(cr < params['cr_high'])):
                        count_rate.append(float(cr))    # per second
                        timeinsec.append(float(sec))
                    else:
                        n_cut += 1
                        continue

                    # put it all together
                    if first_one == 0:
                        parab_phase = parab_phase_sub
                        baseline = baseline_sub
                        timestamp = times
                        first_one = 1
                    else:
                        parab_phase = np.concatenate((parab_phase,parab_phase_sub))
                        baseline = np.concatenate((baseline,baseline_sub))
                        timestamp = np.concatenate((timestamp,times))


                if (n_cut == len(data.getNode(pstring))):
                    failure(row, rarray, 1)
                    continue

                timestamp = np.array(timestamp,dtype=float)
                parab_phase = np.array(parab_phase, dtype=float)
                baseline = np.array(baseline, dtype=float)
                timeinsec = np.array(timeinsec, dtype=float)
                count_rate = np.array(count_rate, dtype=float)

                # subtract out the baseline
                parab_phase -= baseline  

                # cut on no data
                rangex = max(parab_phase)-min(parab_phase)
                if rangex == 0.0:
                    failure(row, rarray, 1)
                    continue


                ## Histogram of the phase amplitudes for each photon for pixel i,j.  0 phase sits at 0 degrees, and a
                ## pulse moves along the loop in the "negative" direction.  So blue produces a bigger response, = more negative value
 
                n_inbin, phasebins = np.histogram(parab_phase, rangex, range = (min(parab_phase),max(parab_phase)))
                binwidth = phasebins[1]-phasebins[0]                    # 1.0
                phasebins += binwidth/2.0
                phasebins = phasebins[0:len(phasebins)-1]

                totalcts = np.sum(n_inbin)
        
                # Cut on low number of total counts
                if totalcts < params['low_counts_val'] * 100.:
                    failure(row, rarray, 1)
                    continue
    
                # Smooth
                window = 'hanning'
                windowsize = int(binwidth * params['bin_smooth'])    
                try:
                    parab_smooth = smooth.smooth(n_inbin, windowsize, window)
                    smoothed_data = np.array(parab_smooth, dtype=float)
                except:
                    failure(row, rarray, 1)
                    continue 

                # Find extreema:
                coarse_data_len = np.floor(len(smoothed_data)/params['big_step'])
                coarse_data = np.zeros(coarse_data_len)
                coarse_x = np.zeros(coarse_data_len)
        
                for l in range(len(coarse_data)):
                    coarse_data[l] = smoothed_data[l*params['big_step']]
                    coarse_x[l] = np.array(phasebins)[l*params['big_step']]

                start_ind = (coarse_data > 2.).nonzero()[0][0]
                coarse_data = coarse_data[start_ind:]
                coarse_x = coarse_x[start_ind:]
        
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
                    
                    if ((grad > 0) & (gradients[count] < 0) & (grad != gradients[count]) & (coarse_data[count] > params['low_counts_val'])):
                        maxima_num += 1
                        max_locations.append(coarse_x[count])
                        max_vals.append(coarse_data[count])
                    if ((grad < 0) & (gradients[count] > 0) & (grad != gradients[count])):
                        minima_num += 1
                        min_locations.append(coarse_x[count])
                        min_vals.append(coarse_data[count])

                if (maxima_num == 2) & (minima_num == 1):
                # set right side minimum
                    minima_num += 1
                    min_vals.append(np.min(coarse_data[np.where(coarse_x > max_locations[-1])]))
                    min_locations.append(coarse_x[np.where(coarse_data == min_vals[-1])][0])

                # Case too many peaks
                if (maxima_num < 2) | (maxima_num > 3) :
                    failure(row, rarray, 2)
                    continue
                if (minima_num > 3):
                    failure(row, rarray, 2)
                    continue
                # Case peak cut off
                if (min_vals[-1] == 0.):
                    failure(row, rarray, 2)
                    continue


                # Find peaks: Tries fitting 2 Gaussians

                ind_right = np.max(np.where(phasebins < min_locations[-1]))
                ind_left = np.min(np.where(phasebins > coarse_x[0]))
                ind_leftg = ind_left

                amplitude1 = max_vals[0]
                amplitude2 = max_vals[1]
 
                x_offset1 = max_locations[0]
                x_offset2 = max_locations[1]

                fwhm = 100.                                   # Could do better...
                sigma1 = fwhm/params['fwhm2sig']
                sigma2 = fwhm/params['fwhm2sig']
        
                paramsin=[sigma1, x_offset1, amplitude1, sigma2, x_offset2, amplitude2]  # First guess at fit params
                errs = np.sqrt(n_inbin)                         # Poisson counts 
                errs[np.where(errs == 0.)] = 1.
                quiet=True

                parinfo = [ {'n':0,'value':paramsin[0],'limits':[fwhm/10., fwhm*2.],             'limited':[True,True],'fixed':False,'parname':"Sigma Blue",'error':0},
                            {'n':1,'value':paramsin[1],'limits':[x_offset1-fwhm, x_offset1+fwhm],'limited':[True,True],'fixed':False,'parname':"x offset blue",'error':0},
                            {'n':2,'value':paramsin[2],'limits':[0., 2.*amplitude1],            'limited':[True,True],'fixed':False,'parname':"Amplitude blue",'error':0},
                            {'n':3,'value':paramsin[3],'limits':[fwhm/10., fwhm*2.],            'limited':[True,True],'fixed':False,'parname':"Sigma red",'error':0},
                            {'n':4,'value':paramsin[4],'limits':[x_offset2-fwhm, x_offset2+fwhm],'limited':[True,True],'fixed':False,'parname':"x offset red",'error':0},
                            {'n':5,'value':paramsin[5],'limits':[0., 2.*amplitude2],            'limited':[True,True],'fixed':False,'parname':"Amplitude red",'error':0}]

                fa = {'x':phasebins[ind_left:ind_right],'y':n_inbin[ind_left:ind_right],'err':errs[ind_left:ind_right]}
                m = mpfit.mpfit(twogaussian, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
            
                mpp = m.params                                #The fit params
                mpperr = m.perror
                chi2gauss = m.fnorm
                redchi2gauss = chi2gauss/len(phasebins[ind_left:ind_right])

                # cut on chi^2
                if (redchi2gauss > params['chi2_cutoff']):
                    failure(row, rarray, 2)
                    continue
        
                for r,p in enumerate(mpp):
                    parinfo[r]['value'] = p
                    if r==0: sigma1 = p
                    if r==1: x_offset1 = p
                    if r==2: amplitude1 = p
                    if r==3: sigma2 = p
                    if r==4: x_offset2 = p
                    if r==5: amplitude2 = p

                # Blobby, got through but cutting now...
                if (np.abs(x_offset2 - min_locations[-1]) < sigma2) | (np.abs(x_offset1 - x_offset2) < 2 * sigma2) | (x_offset2 > min_locations[-1]):
                    failure(row, rarray, 2)
                    continue

                gaussfit = amplitude1 * np.exp( -(pow((phasebins-x_offset1),2) / (2. * pow(sigma1,2)))) + amplitude2 * np.exp( -(pow((phasebins-x_offset2),2) / (2. * pow(sigma2,2)))) 

                # Now fit parabola to get wavelength <-> phase amplitude correspondance
                wavelengths = [params['bluelambda'],params['redlambda']]
                energies = [params['h'] * params['c'] / (x * params['ang2m']) for x in wavelengths]
                energies.append(0.)
                phaseamps = [x_offset1,x_offset2]
                phaseamps.append(0.)

                # Initial guess
                x_off = 0.
                y_off = 0.
                amp = energies[0] / (pow(phaseamps[0],2))

                paramsin=[x_off, amp, y_off]  # First guess at fit params
                errs = np.array([1.,1.,1.])                        
                quiet=True

                parinfo = [ {'n':0,'value':paramsin[0],'limits':[0,0],'limited':[False,False],'fixed':False,'parname':"parabola x offset",'error':0},
                           {'n':1,'value':paramsin[1],'limits':[0,0],'limited':[False,False],'fixed':False,'parname':"parabola amplitude",'error':0},
                           {'n':2,'value':paramsin[2],'limits':[0,0],'limited':[False,False],'fixed':False,'parname':"parabola y offset",'error':0}]

                fa = {'x':phaseamps,'y':energies,'err':errs}
                m = mpfit.mpfit(parabola, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
            
                mpp = m.params                                #The fit params
                mpperr = m.perror
                chi2 = m.fnorm
                redchi2 = chi2/len(phaseamps)
        
                for r,p in enumerate(mpp):
                    parinfo[r]['value'] = p
                    if r==0: x_off = p
                    if r==1: amp = p
                    if r==2: y_off = p

                fitfunc = amp * (pow((phaseamps - x_off), 2 )) + y_off

                # Wavelength/E spectrum
                e_fromphase = amp * (pow((phasebins - x_off), 2 )) + y_off
                lambda_fromphase = (params['h'] * params['c'] / params['ang2m']) / e_fromphase
            
                # Get sigma and wavelength cutoffs:
                
                blue_peak = amp * (pow((x_offset1 - x_off), 2 )) + y_off
                red_peak = amp * (pow((x_offset2 - x_off), 2 )) + y_off

                ind_blue = (np.where(e_fromphase < blue_peak))[0][0]
                ind_red = (np.where(e_fromphase < red_peak))[0][0]

                blue_amp = np.mean(n_inbin[ind_blue-10:ind_blue+10])
                red_amp = np.mean(n_inbin[ind_red-10:ind_red+10])
   
                blue_sigma = np.abs(amp * (pow((x_offset1 - x_off), 2 ) - pow(((x_offset1 + sigma1) - x_off), 2 )))
                red_sigma = np.abs(amp * (pow((x_offset2 - x_off), 2 ) - pow(((x_offset2 + sigma2) - x_off), 2 )))

                try:
                    ind_right = (np.where(e_fromphase < blue_peak + 2.*blue_sigma))[0][0]
                except:
                    ind_right = 0
                try:
                    ind_left = (np.where(e_fromphase < red_peak - 1.5*red_sigma))[0][0]
                except:
                    ind_left = ind_red + 50
        
                paramsin = [blue_sigma, blue_peak, blue_amp, red_sigma, red_peak, red_amp]  # First guess at fit params
                errs = np.sqrt(n_inbin)                         # Poisson counts 
                errs[np.where(errs == 0.)] = 1.
                quiet=True

                parinfo = [ {'n':0,'value':paramsin[0],'limits':[blue_sigma/2.,blue_sigma*1.5],              'limited':[True,True],'fixed':False,'parname':"Sigma Blue",'error':0},
                           {'n':1,'value':paramsin[1],'limits':[blue_peak-blue_sigma, blue_peak+blue_sigma],'limited':[True,True],'fixed':False,'parname':"x offset blue",'error':0},
                           {'n':2,'value':paramsin[2],'limits':[blue_amp/2., blue_amp*2.],                  'limited':[True,True],'fixed':False,'parname':"Amplitude blue",'error':0},
                           {'n':3,'value':paramsin[3],'limits':[red_sigma/2.,red_sigma*1.5],                'limited':[True,True],'fixed':False,'parname':"Sigma red",'error':0},
                           {'n':4,'value':paramsin[4],'limits':[red_peak-red_sigma, red_peak+red_sigma],    'limited':[True,True],'fixed':False,'parname':"x offset red",'error':0},
                           {'n':5,'value':paramsin[5],'limits':[red_amp/2., red_amp*2.],                  'limited':[True,True],'fixed':False,'parname':"Amplitude red",'error':0}]

                fa = {'x':e_fromphase[ind_left:ind_right],'y':n_inbin[ind_left:ind_right],'err':errs[ind_left:ind_right]}
                m = mpfit.mpfit(twogaussian, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
            
                mpp = m.params                                #The fit params
                mpperr = m.perror
                chi2e = m.fnorm
                redchi2e = chi2e/len(e_fromphase[ind_right:ind_left])
        
                for r,p in enumerate(mpp):
                    parinfo[r]['value'] = p
                    if r==0: blue_sigma = p
                    if r==1: blue_peak = p
                    if r==2: blue_amp = p
                    if r==3: red_sigma = p
                    if r==4: red_peak = p
                    if r==5: red_amp = p

                gaussEfit = blue_amp * np.exp( -(pow((e_fromphase-blue_peak),2) / (2. * pow(blue_sigma,2)))) + amplitude2 * np.exp( -(pow((e_fromphase-red_peak),2) / (2. * pow(red_sigma,2)))) 

                # Wavelength range:
                lambda_start = lambda_fromphase[0]

                gradients = np.diff(parab_smooth)

                n_max = 0
                n_min = 0
                loc_max = []
                loc_min = []
                lambda_min = []
                val_max = []
                val_min = []
                count = 0
                for grad in gradients[:-1]:
                    count += 1

                    if ((grad > 0) & (gradients[count] < 0) & (grad != gradients[count]) & (parab_smooth[count] > params['low_counts_val'])):
                        n_max += 1
                        loc_max.append(e_fromphase[count])
                        val_max.append(parab_smooth[count])
                    if ((grad < 0) & (gradients[count] > 0) & (grad != gradients[count])):
                        n_min += 1
                        loc_min.append(e_fromphase[count])
                        lambda_min.append(lambda_fromphase[count])
                        val_min.append(parab_smooth[count])

                lambda_stop = lambda_min[-1]


                # Success - write out solution
                if row['wave_flag'] == 0:
                    row['polyfit'] = np.array([x_off, y_off, amp])
                    row['sigma'] = blue_sigma                                 # sigma of blue peak in eV
                    row['solnrange'] = np.array([lambda_start,lambda_stop])   # Angstroms
                    driftrow['pixelrow'] = i
                    driftrow['pixelcol'] = j
                    driftrow['gaussparams'] = np.array([x_offset1, amplitude1, sigma1, x_offset2, amplitude2, sigma2])
                    driftrow.append()
                    resest = np.abs(blue_peak) / (params['fwhm2sig'] * blue_sigma)
                    rarray.append(resest)

                    # fits plot
                    if (plotcounter % n_plots_per_page == 0):
                        fig6 = plt.figure(figsize=(8.25, 10), dpi=100)

                    plt.subplot(5,4,plotcounter%20+1)
                    plt.plot(phasebins,n_inbin, label='raw')
                    titlestring = str(i)+', '+str(j)+' '+pstring
                    plt.title(titlestring)
                    plt.xlim(phasebins[ind_leftg]-50.,np.min([-50.,np.max(phasebins)]))
                    plt.ylim(0.,max_vals[-1]+60.)
                    plt.plot(phasebins, gaussfit, 'r', label='Gaussian fit1')
        
                    if (((plotcounter +1) % n_plots_per_page == 0) or ((i+1 == n_rows) and (j+1 == n_cols))):
                        pp.savefig(fig6)
            
                    plotcounter+=1

                else:
                    failure(row, rarray, 2)
                    continue


            # insert row; writes to the table I/O buffer
            row.append()

        # flush the table's I/O buffer to write the data to disk
        caltable.flush()
        drifttable.flush()

        # close the file, flush all remaining buffers
        h5out.close()
        h5outdrift.close()

        data.close()

        pp.close()
        print plotcounter, ' good pixels'

        # Plots

        plotArray( xarray, yarray, rarray, colormap=mpl.cm.gnuplot2, showMe=False,
              plotFileName=outdir+datedir+params['figdir']+outfile.split('.')[0]+'_arrayPlot.png', plotTitle='Energy Resolution')

        n_res, resbins = np.histogram(rarray, 80, range = (1,12))
        binwidth = (resbins[1]-resbins[0])/2.0
        resbins += binwidth
        fig3=plt.figure()
        ax2=fig3.add_subplot(111)
        ax2.plot(resbins[:-1],n_res)
        plt.xlabel('Energy Resolution')
        plt.savefig(outdir+datedir+params['figdir']+outfile.split('.')[0]+"_R_Estimates.png")
        plt.clf()
        

    print 'Finish Time: ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

   

if __name__ == '__main__':
    
    paramFile = sys.argv[1]
    wavelengthCal(paramFile)






