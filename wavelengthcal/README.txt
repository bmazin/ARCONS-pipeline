Danica Marsden                                                                           November 1, 2012

README file for the wavelength calibration module:

Wavelength calibration reads in a wavelength calibration file ('cal_...h5'; blue, red and IR laser 
light) that has a beam map solution, and makes a histogram of photon phase amplitudes for 
each pixel. Bad pixels are flagged.  Gaussians are fit to the histograms.

** The wavelength file must be contained in a file list, which can be
   generated with list_calfiles.sh and produces calfile_list.txt.  The
   output directories should match.

The peaks are fit with a polynomial and that gives a phase amplitude <-> wavelength 
correspondence.


For the Lick 2012 run:
blue = 400.94 nm
red  = 656.04 nm
IR   = 974.49 nm


Dependencies (in the ARCONS-pipeline utils directory):
- readDict.py
- mpfit.py
- smooth.py
- fitFunctions.py


To run inside another script:

paramFile = sys.argv[1]
wavelengthCal(paramFile)

where paramFile has the same parameters as e.g. waveCal.dict from the ARCONS-pipeline
params directory


The output file, 'calsol_...h5' is a table where each row is populated with a WaveCalSoln:

class WaveCalSoln(IsDescription):
    roach = tables.UInt16Col()              # ROACH board number
    pixelnum = tables.UInt16Col()           # pixel number on the roach
    pixelrow = tables.UInt16Col()           # physical x location - from beam map
    pixelcol = tables.UInt16Col()           # physical y location 
    polyfit = tables.Float64Col(3)          # polynomial to convert from phase amplitude to wavelength,
                                                            #    double precision
    sigma = tables.Float64Col()             # 1 sigma (Gaussian width) in eV
    solnrange = tables.Float32Col(2)      # start and stop wavelengths for
    	      				                      #  the fit in Angstroms
    wave_flag = tables.UInt16Col()          # flag to indicate if pixel is good (0), dead (1) or failed
                                                             #    during wave cal fitting (2)
  
where the table is created in the following way:

h5out = tables.openFile(outfile, mode='w')
root = h5out.root
calgroup = h5out.createGroup(root, 'wavecal', 'Table of calibration parameters for each pixel')
caltable = h5out.createTable(calgroup, 'calsoln', WaveCalSoln, title='Wavelength Cal Table')
row = caltable.row
for i in range(n_rows):
    for j in range(n_cols):
        row['pixelrow'] = i 
        row['pixelcol'] = j
 
          ... etc.
        row.append()
caltable.flush()
h5out.close()


Outputs:

The polynomial fit parameters for each pixel, are for a parabola:
3 coefficients -  x_offset, y_offset, amplitude such that

Energy in eV = amplitude * (x - x_offset)^2 + y_offset

and then wavelength in angstroms = ( h [eV s] * c [m/s] / 1.e-10 [ang/m] )   / Energy in eV

The enegry error is in eV, and corresponds to the Gaussian sigma when fitting the 
    counts histogram at the blue peak.
The start and stop wavelengths indicate the range over which the
    solution can be trusted, in Angstroms.
If a pixel is good, the flag is set to 0, otherwise it is set to another integer (1 for a dead or defunct 
    pixel, 2 for a bad fit based on chi^2) and the polyfit entries are set to -1. 
