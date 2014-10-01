Danica Marsden                                                                           November 1, 2012
Alex Walter                                                                              August 4, 2014

README file for the wavelength calibration module:

SUMMARY:
    Wavelength calibration reads in a wavelength calibration file ('cal_...h5'; blue, red and IR laser 
    light) that has a beam map solution, and makes a histogram of photon phase amplitudes for 
    each pixel. Bad pixels are flagged.  Gaussians are fit to the histograms. If the 3 laser gaussian fit fails 
    it tries to just fit the blue and red laser peaks instead. The peaks are fit with a polynomial and that 
    gives a phase amplitude <-> wavelength correspondence.

USAGE:
    The waveCal class is initiated with a FileName object and the waveCal.dict parameter dictionary: 
        waveCalObject=waveCal(calFN,params,save_pdf=True,verbose=False,debug=False)
    To find the wavelength solutions you call findWaveLengthSoln(). If save_pdf=True this saves a pdf of the phase histograms in the output figs folder.
        waveCalObject.findWaveLengthSoln()
    To write the wavecal solution to an H5 file and the _drift file you call:
        waveCalObject.write_waveCal()
        waveCalObject.write_waveCal_drift()
    Once the H5 files exist you can make further diagnostic plots with the waveCal_diagnostic class:
        diag_obj=waveCal_diagnostic(calFN,params,save=True)
        diag_obj.make_R_array()
        diag_obj.plot_R_array()
        diag_obj.plot_nlaser_array()
        diag_obj.plot_R_hist()

    Often we want to find wavecal solutions on a number of cal files. To do this you just loop through a list of them. 
    A list of cal FileNames can be generated with getCalFileNames() using the waveCal.dict parameter dictionary:
        calFNs, params = getCalFileNames(paramFile)

INPUT:
    calFN - FileName object of cal file to be analyzed
    params - dictionary of parameters (see /params/waveCal.dict and util.readDict)
    save_pdf - if True, saves a pdf of phase histograms for each non-dead pixel
    verbose - set True for runtime comments
    debug - set True to raise errors (most pixel failures). Automatically plots and pauses at each pixel solution
    
OUTPUT:
    - wavecal solution H5 file (see below)
    - wavecal _drift H5 file (see below)
    - pdf of phase histograms

DEPENDENCIES:
    - util/readDict.py
    - util/mpfit.py
    - util/smooth.py
    - util/fitFunctions.py
    - util/ObsFile.py
    - util/FileName.py
    - hotpix/hotPixels.py


WAVECAL SOLTUION FILE:
    The output file, 'calsol_...h5' is a table where each row is populated with a WaveCalSoln:

    WaveCalSoln_Description = {
        "roach"     : UInt16Col(),      # ROACH board number
        "pixelnum"  : UInt16Col(),      # pixel number on the roach
        "pixelrow"  : UInt16Col(),      # physical x location - from beam map
        "pixelcol"  : UInt16Col(),      # physical y location 
        "polyfit"   : Float64Col(3),    # polynomial to convert from phase amplitude to wavelength float 64 precision
        "sigma"     : Float64Col(),     # 1 sigma (Gaussian width) in eV, for blue peak
        "solnrange" : Float32Col(2),    # start and stop wavelengths for the fit in Angstroms
        "wave_flag" : UInt16Col()}      # flag to indicate if pixel is good (0), unallocated (1), dead (2), or failed during wave cal fitting (2+)   

    The 3 polynomial fit parameters for each pixel, are for a parabola:
        polyfit = [constant_term, linear_term, quadratic_term]
        Energy (eV) = constant_term + x*linear_term + x^2*quadratic_term
    where x is the phase in ADC/DAC units.

    The enegry error is in eV, and corresponds to the Gaussian sigma when fitting the counts histogram at the blue peak.

    The start and stop wavelengths indicate the range over which the solution can be trusted, in Angstroms.

    wave_flags:
    0 - Success!
    1 - Pixel not in beammap
    2 - Dead pixel or count rate less than 'min_count_rate'
    3 - Unable to find blue laser peak
    4 - Fit failed on blue laser peak
    5 - Blue peak fit has chi^2 larger than 'max_chi2_blue'
    6 - Unable to guess parameters for blue/red/IR/noise fit
    ====== If 3 laser peak fails, try 2 laser peak =====
    7 - Unable to find blue/red laser peaks
    8 - Fit failed on blue/red laser peaks
    9 - Fit hit parameter limits
    10 - 
    11 -
    12 - Fit has chi^2 larger than 'max_chi2_all'
    13 - Parabola fit failed


WAVECAL DRIFT FILE:
    The output file, 'calsol_..._drift.h5' is a table where successful pixels have their fit parameters saved. If no pixels have a successful fit then no drift file is made

    DriftObj_Description = {
        "pixelrow"      : UInt16Col(),                  # physical x location - from beam map
        "pixelcol"      : UInt16Col(),                  # physical y location 
        "gaussparams"   : Float64Col(num_params),        # parameters used to fit data
        "perrors"       : Float64Col(num_params)}        # the errors on the fits
        
    gaussparams are the parameters for the phase histogram gaussians and noise tail:
    p[0] = sigma1
    p[1] = x_offset1
    p[2] = amplitude1
    p[3] = sigma2
    p[4] = x_offset2
    p[5] = amplitude2
    p[6] = sigma3
    p[7] = x_offset3
    p[8] = amplitude3
    p[9] = scale_factor4
    p[10] = x_offset4
    p[11] = amplitude4
    gauss1 = p[2] * np.exp( - (pow(( x - p[1]),2) / ( 2. * pow(p[0],2))))
    gauss2 = p[5] * np.exp( - (pow(( x - p[4]),2) / ( 2. * pow(p[3],2))))
    gauss3 = p[8] * np.exp( - (pow(( x - p[7]),2) / ( 2. * pow(p[6],2))))
    power4 = p[11] * np.maximum((x - p[10]),0)**p[9]
    model = gauss1 + gauss2 + gauss3 + power4

