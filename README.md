ARCONS-pipeline
===============

This repository contains the code to analyze data from the ARCONS array.  

***

Required external software components:
---------------------

Enthought Python Distribution (EPD) 7.3 (http://www.enthought.com/products/epd.php)
 
PyEphem (http://rhodesmill.org/pyephem/) 

(you can check if you have them with help('modules') within the (i)python interpreter)

PyTables (see http://www.tumblr.com/tagged/pytables and instructions therein for Mac) 

***

Recommended external software components:
---------------------

Aptana Studio 3 (http://www.aptana.com/products/studio3/download)

***

General purpose:
---------------------

/headers contains standard definitions 

/params contains files that provide inputs to the pipeline 

/utils contains commonly used functions 

/examples contains simple examples to show how to use the software


Pipeline components:
---------------------

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

***

This document uses the markdown syntax, see http://daringfireball.net/projects/markdown/
