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
    
def failure(row, rarray, larray, flagnum):
    row['wave_flag'] = flagnum
    row['polyfit'] = np.array([-1.,-1., -1.])
    row['sigma'] = -1.
    row['solnrange'] = np.array([-1.,-1.])
    row.append()
    rarray.append(0.)
    larray.append(0.)

def fitTwo(peaks, peak_locations, xarr, yarr):

    amplitude1 = peaks[0]
    amplitude2 = peaks[1]
 
    x_offset1 = peak_locations[0]
    x_offset2 = peak_locations[1]

    fwhm = 100.                                   # Could do better...
    fwhm2sig = 2.355
    sigma1 = fwhm/fwhm2sig
    sigma2 = fwhm/fwhm2sig
        
    params2=[sigma1, x_offset1, amplitude1, sigma2, x_offset2, amplitude2]  # First guess at fit params
    errs = np.sqrt(yarr)                         # Poisson counts 
    errs[np.where(errs == 0.)] = 1.
    quiet=True

    parinfo = [ {'n':0,'value':params2[0],'limits':[fwhm/10., fwhm*2.],            'limited':[True,True],'fixed':False,'parname':"Sigma Blue",'error':0},
               {'n':1,'value':params2[1],'limits':[x_offset1-fwhm, x_offset1+fwhm],'limited':[True,True],'fixed':False,'parname':"x offset blue",'error':0},
               {'n':2,'value':params2[2],'limits':[0., 2.*amplitude1],            'limited':[True,True],'fixed':False,'parname':"Amplitude blue",'error':0},
               {'n':3,'value':params2[3],'limits':[fwhm/10., fwhm*2.],            'limited':[True,True],'fixed':False,'parname':"Sigma red",'error':0},
               {'n':4,'value':params2[4],'limits':[x_offset2-fwhm, x_offset2+fwhm],'limited':[True,True],'fixed':False,'parname':"x offset red",'error':0},
               {'n':5,'value':params2[5],'limits':[0., 2.*amplitude2],            'limited':[True,True],'fixed':False,'parname':"Amplitude red",'error':0}]

    fa = {'x':xarr,'y':yarr,'err':errs}

    m = mpfit.mpfit(twogaussian, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
            
    mpp = m.params                                #The fit params
    mpperr = m.perror
    chi2gauss = m.fnorm
    redchi2gauss2 = chi2gauss/len(xarr)
        
    for k,p in enumerate(mpp):
        parinfo[k]['value'] = p
        #print parinfo[k]['parname'],p," +/- ",mpperr[j]
        if k==0: sigma1 = p
        if k==1: x_offset1 = p
        if k==2: amplitude1 = p
        if k==3: sigma2 = p
        if k==4: x_offset2 = p
        if k==5: amplitude2 = p
    params2=[sigma1, x_offset1, amplitude1, sigma2, x_offset2, amplitude2]
    
    return params2, redchi2gauss2
        

def fitThree(peaks, peak_locations, xarr, yarr):

    amplitude1 = peaks[0]
    amplitude2 = peaks[1]
    amplitude3 = peaks[2]
 
    x_offset1 = peak_locations[0]
    x_offset2 = peak_locations[1]
    x_offset3 = peak_locations[2]

    fwhm = 100.
    fwhm2sig = 2.355
    sigma1 = fwhm/fwhm2sig
    sigma2 = fwhm/fwhm2sig
    sigma3 = fwhm/fwhm2sig
        
    params3=[sigma1, x_offset1, amplitude1, sigma2, x_offset2, amplitude2, sigma3, x_offset3, amplitude3]  
    errs = np.sqrt(yarr)                         
    errs[np.where(errs == 0.)] = 1.
    quiet=True

    parinfo = [ {'n':0,'value':params3[0],'limits':[fwhm/10., fwhm*2.],            'limited':[True,True],'fixed':False,'parname':"Sigma Blue",'error':0},
               {'n':1,'value':params3[1],'limits':[x_offset1-fwhm, x_offset1+fwhm],'limited':[True,True],'fixed':False,'parname':"x offset blue",'error':0},
               {'n':2,'value':params3[2],'limits':[0., 2.*amplitude1],            'limited':[True,True],'fixed':False,'parname':"Amplitude blue",'error':0},
               {'n':3,'value':params3[3],'limits':[fwhm/10., fwhm*2.],            'limited':[True,True],'fixed':False,'parname':"Sigma red",'error':0},
               {'n':4,'value':params3[4],'limits':[x_offset2-fwhm, x_offset2+fwhm],'limited':[True,True],'fixed':False,'parname':"x offset red",'error':0},
               {'n':5,'value':params3[5],'limits':[0., 2.*amplitude2],            'limited':[True,True],'fixed':False,'parname':"Amplitude red",'error':0},
               {'n':6,'value':params3[6],'limits':[fwhm/10., fwhm*2.],            'limited':[True,True],'fixed':False,'parname':"Sigma ir",'error':0},
               {'n':7,'value':params3[7],'limits':[x_offset3-fwhm, x_offset3+fwhm],'limited':[True,True],'fixed':False,'parname':"x offset ir",'error':0},
               {'n':8,'value':params3[8],'limits':[0., 2.*amplitude3],            'limited':[True,True],'fixed':False,'parname':"Amplitude ir",'error':0}]

    fa = {'x':xarr,'y':yarr,'err':errs}

    m = mpfit.mpfit(threegaussian, functkw=fa, parinfo=parinfo, maxiter=1000, quiet=quiet)
            
    mpp = m.params                                
    mpperr = m.perror
    chi2gauss = m.fnorm
    redchi2gauss3 = chi2gauss/len(xarr)
        
    for k,p in enumerate(mpp):
        parinfo[k]['value'] = p
        if k==0: sigma1 = p
        if k==1: x_offset1 = p
        if k==2: amplitude1 = p
        if k==3: sigma2 = p
        if k==4: x_offset2 = p
        if k==5: amplitude2 = p    
        if k==6: sigma3 = p
        if k==7: x_offset3 = p
        if k==8: amplitude3 = p   
    params3=[sigma1, x_offset1, amplitude1, sigma2, x_offset2, amplitude2, sigma3, x_offset3, amplitude3]

    return params3, redchi2gauss3

    
def wavelengthCal(paramFile): 


    params = readDict.readDict(paramFile)
    params.readFromFile(paramFile)

    infile = params['calfileslist']
    files2run = params['files2run']
    outdir = params['outdir']

    cal_files = []

    file_list = open(outdir+infile, 'r')

    if (('cal_' in files2run) | ('obs_' in files2run)):
        cal_files.append(files2run)
    else:
        if '2012' in files2run:
            for line in file_list:
                if files2run in line.split('cal')[0]:
                    cal_files.append(line.strip())
        else:
            for line in file_list:
                cal_files.append(line.strip())
    

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
        larray = []
        roacharr = []

        # Open outfile to write parameters to
        try:
            outfile = 'calsol' + cal_files[k].split('cal')[1]
        except:
            outfile = 'calsol' + cal_files[k].split('obs')[1]
            
        try:
            h5out = tables.openFile(outdir+datedir+outfile, mode='w')
        except:
            pass
        root = h5out.root

        # Create a group that branches from the root node, save cal data table in this group
        calgroup = h5out.createGroup(root, 'wavecal', 'Table of calibration parameters for each pixel')

        # Create a table instance under group calgroup with node name "calsoln"
        # the WaveCalSoln class declared before is the description parameter to define the columns of the table
        caltable = h5out.createTable(calgroup, 'calsoln', WaveCalSoln, title='Wavelength Cal Table')

        # Pointer to the table row instance:
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


        # Get the beam map, an array of strings that translates physical location [row_i,col_j] to the roach/pixel in the readout
        beammap = data.root.beammap.beamimage.read()  


        # Plot of fits for each "good" pixel
        try:
            os.mkdir(outdir+datedir+params['figdir'])
        except:
            pass
        
        pp = PdfPages(outdir+datedir+params['figdir']+outfile.split('.')[0]+"_fits.pdf")
        mpl.rcParams['font.size'] = 4
        n_plots_per_page = 20
        plotcounter = 0


        # Loop through each pixel and get its solution, populate output WaveCalSoln and append the row.
        for i in range(n_rows):
            for j in range(n_cols):

                pstring = beammap[i][j]
                pixelNode = pstring.split('t')[0]
    
                xarray.append(j)
                yarray.append(i)
                roach = int(pixelNode.split('/')[1][1:])
                pixel = int(pixelNode.split('/')[2][1:])
                roacharr.append(roach)

                row['pixelrow'] = i
                row['pixelcol'] = j
                row['roach'] = roach
                row['pixelnum'] = pixel
                row['wave_flag'] = 0                      # 0 until flagged (1 for dead, 2 for bad data/fit)

                # The non-sorted list of all photon packets belonging to pixel i,j
                photondata = np.concatenate(data.root._f_getChild(pstring)[:])

                # Flag dead pixels
                if len(photondata)==0:
                    failure(row, rarray, larray, 1)
                    continue

                # For each second, get the packets
                timeinsec = []
                count_rate = []
                first_one = 0
                for sec in range(len(data.getNode(pstring))):
            
                    packets = data.getNode(pstring)[sec]
            
                    microtimes = packets & params['timeMask']
                    parab_phase_sub = packets >> params['n_bits_parab'] & params['pulseMask']
                    baseline_sub = packets >> params['n_bits_base'] & params['pulseMask']

                    times = float(sec) + microtimes*(1.e-6)

                    # Cut on count rate - data sliced out where pixels go hot
                    cr = len(packets)
                    count_rate.append(float(cr))    # per second
                    timeinsec.append(float(sec))

                    # Put it all together
                    if first_one == 0:
                        parab_phase = parab_phase_sub
                        baseline = baseline_sub
                        timestamp = times
                        first_one = 1
                    else:
                        parab_phase = np.concatenate((parab_phase,parab_phase_sub))
                        baseline = np.concatenate((baseline,baseline_sub))
                        timestamp = np.concatenate((timestamp,times))

                timestamp = np.array(timestamp,dtype=float)
                parab_phase = np.array(parab_phase, dtype=float)
                baseline = np.array(baseline, dtype=float)
                timeinsec = np.array(timeinsec, dtype=float)
                count_rate = np.array(count_rate, dtype=float)

                # Subtract out the baseline
                parab_phase -= baseline  

                # Chop out high count rate bits:
                rate_high = np.where(count_rate > np.mean(count_rate) + params['cr_high'])
                rate_ok = np.where(count_rate <= np.mean(count_rate) + params['cr_high'])
                if len(rate_high[0]) != 0:
                    goodsec = timeinsec[rate_ok[0]]
                    goodind = []
                    for q in range(len(timestamp)):
                        if (np.floor(timestamp[q]) in goodsec):
                            goodind.append(q)
                    parab_phase = parab_phase[goodind]

                # Cut off any hits above zero (~nonsensible)
                parab_phase = parab_phase[np.where(parab_phase < 0.)[0]]

                print i, j, parab_phase
                # Cut on no data
                rangex = max(parab_phase)-min(parab_phase)
                if rangex == 0.0:
                    failure(row, rarray, larray, 1)
                    continue


                # Histogram of the phase amplitudes for each photon for pixel i,j.  0 phase sits at 0 degrees, and a
                # pulse moves along the loop in the "negative" direction.  So blue produces a bigger response, = more negative value
 
                n_inbin, phasebins = np.histogram(parab_phase, rangex, range = (min(parab_phase),max(parab_phase)))
                binwidth = phasebins[1]-phasebins[0]                    # 1.0
                phasebins += binwidth/2.0
                phasebins = phasebins[0:len(phasebins)-1]

                totalcts = np.sum(n_inbin)
        
                # Cut on low number of total counts
                if totalcts < params['low_counts_val'] * 100.:
                    failure(row, rarray, larray, 1)
                    continue
    
                # Smooth
                window = 'hanning'
                windowsize = int(binwidth * params['bin_smooth'])    
                try:
                    parab_smooth = smooth.smooth(n_inbin, windowsize, window)
                    smoothed_data = np.array(parab_smooth, dtype=float)
                except:
                    failure(row, rarray, larray, 1)
                    continue 

                # Find extreema:
                coarse_data_len = np.floor(len(smoothed_data)/params['big_step'])
                coarse_data = np.zeros(coarse_data_len)
                coarse_x = np.zeros(coarse_data_len)
        
                for l in range(len(coarse_data)):
                    coarse_data[l] = smoothed_data[l*params['big_step']]
                    coarse_x[l] = np.array(phasebins)[l*params['big_step']]

                try:
                    start_ind = (coarse_data > params['low_counts_val']).nonzero()[0][0]
                    end_ind = (coarse_data > params['low_counts_val']).nonzero()[0][-1]
                except:
                    failure(row, rarray, larray, 1)
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
                    
                    if ((grad > 0) & (gradients[count] < 0) & (grad != gradients[count]) & (coarse_data[count] > params['low_counts_val'])):
                        maxima_num += 1
                        max_locations.append(coarse_x[count])
                        max_vals.append(coarse_data[count])
                    if ((grad < 0) & (gradients[count] > 0) & (grad != gradients[count])):
                        minima_num += 1
                        min_locations.append(coarse_x[count])
                        min_vals.append(coarse_data[count])

                # Set right side minimum if more maxima than minima
                if (((maxima_num ==  params['n_lasers']-1) & (minima_num ==  params['n_lasers']-2))
                    | ((maxima_num ==  params['n_lasers']) & (minima_num == params['n_lasers']-1))):
                    minima_num += 1
                    min_vals.append(np.min(coarse_data[np.where(coarse_x > max_locations[-1])]))
                    min_locations.append(coarse_x[np.where(coarse_data == min_vals[-1])][0])

                 # Remove right side maximum if 4 maxima, 3 minima
                if ((maxima_num == params['n_lasers']+1) & (minima_num == params['n_lasers'])):
                    maxima_num -= 1
                    max_vals = max_vals[:params['n_lasers']-1]
                    max_locations = max_locations[:params['n_lasers']-1]          


                # Cut if wrong number of peaks
                if (maxima_num < params['n_lasers']-1) | (maxima_num > params['n_lasers']) :
                    failure(row, rarray, larray, 2)
                    continue
                if (minima_num != maxima_num):
                    failure(row, rarray, larray, 2)
                    continue


                # Peak fitting:

                ind_right = np.max(np.where(phasebins < min_locations[-1]))
                ind_left = np.min(np.where(phasebins > coarse_x[0]))
                ind_leftg = ind_left

                n_in_fit = 0
                # If 3, pass in max_vals, max_locs, phasebins[ind_left:ind_right], n_inbin[ind_left:ind_right] -> fitThree(), return params
                if ((maxima_num == params['n_lasers']) & (len(max_locations) == params['n_lasers'])):
                    n_in_fit = params['n_lasers']
                    gparams, redchi2gauss = fitThree(max_vals, max_locations, phasebins[ind_left:ind_right], n_inbin[ind_left:ind_right]) 
                    sigma1 = gparams[0]
                    x_offset1 = gparams[1]
                    amplitude1 = gparams[2]
                    sigma2 = gparams[3]
                    x_offset2 = gparams[4]
                    amplitude2 = gparams[5]
                    sigma3 = gparams[6]
                    x_offset3 = gparams[7]
                    amplitude3 = gparams[8]
                else:
                    n_in_fit = params['n_lasers']-1
                    gparams, redchi2gauss = fitTwo(max_vals, max_locations, phasebins[ind_left:ind_right], n_inbin[ind_left:ind_right]) 
                    sigma1 = gparams[0]
                    x_offset1 = gparams[1]
                    amplitude1 = gparams[2]
                    sigma2 = gparams[3]
                    x_offset2 = gparams[4]
                    amplitude2 = gparams[5]
            
                if (n_in_fit == params['n_lasers']-1):
                    gaussfit = amplitude1 * np.exp( -(pow((phasebins-x_offset1),2) / (2. * pow(sigma1,2)))) + \
                               amplitude2 * np.exp( -(pow((phasebins-x_offset2),2) / (2. * pow(sigma2,2)))) 
                else:
                    gaussfit = amplitude1 * np.exp( -(pow((phasebins-x_offset1),2) / (2. * pow(sigma1,2)))) + \
                               amplitude2 * np.exp( -(pow((phasebins-x_offset2),2) / (2. * pow(sigma2,2)))) + \
                               amplitude3 * np.exp( -(pow((phasebins-x_offset3),2) / (2. * pow(sigma3,2))))
 

                ## Now fit parabola to get wavelength <-> phase amplitude correspondance

                if (n_in_fit == params['n_lasers']):
                    wavelengths = [params['bluelambda'],params['redlambda'], params['irlambda']]
                    phaseamps = [x_offset1, x_offset2, x_offset3]
                else:
                    wavelengths = [params['bluelambda'],params['redlambda']]
                    phaseamps = [x_offset1, x_offset2]
                phaseamps.append(0.) 
                energies = [params['h'] * params['c'] / (x * params['ang2m']) for x in wavelengths]
                energies.append(0.)

                # Initial guess
                x_off = 0.
                y_off = 0.
                amp = energies[0] / (pow(phaseamps[0],2))

                paramsin=[x_off, amp, y_off]  # First guess at fit params
                errs = np.ones((len(energies)), float)                        
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


                # Final cuts
            
                if (np.abs(gparams[-2] - min_locations[-1]) < gparams[-3]/2.) | (np.abs(gparams[-5] - gparams[-2]) < 2 * gparams[-3]) | (gparams[-2] > min_locations[-1]):
                    failure(row, rarray, larray, 2)
                    continue
                #if (redchi2gauss > params['chi2_cutoff']):              # Cut on chi^2
                #    failure(row, rarray, larray, 2)
                #    continue


                # Wavelength/E spectrum
                
                e_fromphase = amp * (pow((phasebins - x_off), 2 )) + y_off
                lambda_fromphase = (params['h'] * params['c'] / params['ang2m']) / e_fromphase
            
                # Get sigma and wavelength cutoffs:
                
                blue_peak = amp * (pow((x_offset1 - x_off), 2 )) + y_off
                red_peak = amp * (pow((x_offset2 - x_off), 2 )) + y_off

                try:
                    ind_blue = (np.where(e_fromphase < blue_peak))[0][0]
                    ind_red = (np.where(e_fromphase < red_peak))[0][0]
                except:
                    failure(row, rarray, larray, 2)
                    continue

                blue_amp = np.mean(n_inbin[ind_blue-10:ind_blue+10])
                red_amp = np.mean(n_inbin[ind_red-10:ind_red+10])
   
                blue_sigma = np.abs(amp * (pow((x_offset1 - x_off), 2 ) - pow(((x_offset1 + sigma1) - x_off), 2 )))
                red_sigma = np.abs(amp * (pow((x_offset2 - x_off), 2 ) - pow(((x_offset2 + sigma2) - x_off), 2 )))

                params_wave=[blue_sigma, blue_peak, blue_amp, red_sigma, red_peak, red_amp]

                try:
                    ind_right = (np.where(e_fromphase < blue_peak + 2.*blue_sigma))[0][0]
                except:
                    ind_right = 0
                
                if (n_in_fit == params['n_lasers']):
                    ir_peak = amp * (pow((x_offset3 - x_off), 2 )) + y_off
                    ind_ir = (np.where(e_fromphase < ir_peak))[0][0]
                    ir_amp = np.mean(n_inbin[ind_ir-10:ind_ir+10])
                    ir_sigma = np.abs(amp * (pow((x_offset3 - x_off), 2 ) - pow(((x_offset3 + sigma3) - x_off), 2 )))
                    params_wave.append(ir_sigma)
                    params_wave.append(ir_peak)
                    params_wave.append(ir_amp)
                    gaussEfit = blue_amp * np.exp( -(pow((e_fromphase-blue_peak),2) / (2. * pow(blue_sigma,2)))) + \
                                red_amp * np.exp( -(pow((e_fromphase-red_peak),2) / (2. * pow(red_sigma,2)))) + \
                                ir_amp * np.exp( -(pow((e_fromphase-ir_peak),2) / (2. * pow(ir_sigma,2))))
                    try:
                        ind_left = (np.where(e_fromphase < ir_peak - 1.5*ir_sigma))[0][0]
                    except:
                        ind_left = ind_ir + 50
                else:
                    gaussEfit = blue_amp * np.exp( -(pow((e_fromphase-blue_peak),2) / (2. * pow(blue_sigma,2)))) + \
                                red_amp * np.exp( -(pow((e_fromphase-red_peak),2) / (2. * pow(red_sigma,2))))
                    try:
                        ind_left = (np.where(e_fromphase < red_peak - 1.5*red_sigma))[0][0]
                    except:
                        ind_left = ind_red + 50

                
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
                    larray.append(n_in_fit)

                    # fits plot
                    if (plotcounter % n_plots_per_page == 0):
                        fig6 = plt.figure(figsize=(8.25, 10), dpi=100)

                    plt.subplot(5,4,plotcounter%20+1)
                    plt.plot(phasebins,n_inbin, label='raw')
                    titlestring = str(i)+', '+str(j)+' '+pstring
                    plt.title(titlestring)
                    plt.xlim(phasebins[ind_leftg]-50., np.max(phasebins))
                    plt.ylim(0.,np.max(max_vals)+60.)
                    plt.plot(phasebins, gaussfit, 'r', label='Gaussian fit1')
        
                    if (((plotcounter +1) % n_plots_per_page == 0) or ((i+1 == n_rows) and (j+1 == n_cols))):
                        pp.savefig(fig6)
            
                    plotcounter+=1

                else:
                    failure(row, rarray, larray, 2)
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

        plotArray( xarray, yarray, larray, colormap=mpl.cm.gnuplot2, showMe=False,
              plotFileName=outdir+datedir+params['figdir']+outfile.split('.')[0]+'_nlaserPlot.png', plotTitle='Number of Lasers for Fit')

        n_res, resbins = np.histogram(rarray, 80, range = (1,12))
        binwidth = (resbins[1]-resbins[0])/2.0
        resbins += binwidth
        fig3=plt.figure()
        ax2=fig3.add_subplot(111)
        #ax2.plot(resbins[:-1],n_res)
        colormap = mpl.cm.gist_ncar
        plt.gca().set_color_cycle([colormap(i) for i in np.linspace(0, 0.9, params['n_roaches'])])
        roacharr = np.array(roacharr)
        rarray = np.array(rarray)
        labels = []
        for iterator in range(params['n_roaches']):
            roach_pix = np.where(roacharr == iterator)[0]
            n_res, resbins = np.histogram(rarray[roach_pix], 80, range = (1,12))
            ax2.plot(resbins[1:-1],n_res[1:])
            labels.append('roach %i' %(iterator))
        plt.legend(labels, loc='upper left')
        plt.xlim(1.,9.)
        plt.xlabel('Energy Resolution at 400nm')
        plt.savefig(outdir+datedir+params['figdir']+outfile.split('.')[0]+"_R_Estimates.png")
        plt.clf()
        

    print 'Finish Time: ', time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())

   

if __name__ == '__main__':
    
    paramFile = sys.argv[1]
    wavelengthCal(paramFile)






