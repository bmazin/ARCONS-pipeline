This repository contains the code to analyze data from the ARCONS array.

General purpose:
/headers contains standard definitions
/params contains files that provide inputs to the pipeline
/utils contains commonly used functions

Pipeline components:
/cosmic contains a module for cosmic ray cleaning
/wavelengthcal contains a module to do wavelength calibration 
/flatcal contains a module for normalizing the QE as a function of wavelength between pixels
/fluxcal contains a module to calibrating the system compared to a standard star
/astrometry contains a module to determine the RA and DEC of each photon
/skysub contains a module to subtract the sky background
/imagecube contains a module to generate a FITS image cube based on an observations (no timing info)
/legacy contains ARCONS analysis code that predates this repository
/quicklook contains tools for quickly looking at ARCONS h5 files
/beammap contains tools for creating, viewing, and modifying beam maps 

Each directory contains a /test subdirectory, where code to test the module will be stored.