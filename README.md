ARCONS-pipeline
===============

This repository contains the code to analyze data from the ARCONS array.  

***

Required external software components:
---------------------

Enthought Python Distribution (EPD) 7.3 (http://www.enthought.com/products/epd.php)
 
PyEphem (http://rhodesmill.org/pyephem/)

PyGuide (http://www.astro.washington.edu/users/rowen/PyGuide/Manual.html)

(you can check if you have them with help('modules') within the (i)python interpreter)

gittle -- represents the state of the git repository in python (easy_install gittle)

ds9 -- need to download pyds9 and install.  Instructions here:
http://hea-www.harvard.edu/saord/ds9/pyds9/

If you are having troubles with PyTables (which you shouldn't since it is built into EPD), see http://www.tumblr.com/tagged/pytables and instructions therein for Mac.

***
To make nicer looking plots, you need the packages configobj and mpltools

  $ easy_install configobj
  $ git clone git@github.com:tonysyu/mpltools.git
  $ cd mpltools
  $ python setup.py install

***
Another package you will need is pyinterval.  It depends on crlibm

NOTE:  the python package interval is not the same thing.  If this
got installed by mistake, things might not work so well.  A tell-tale
sign that you have the correct thing installed is that this line 
works in ipython:

from interval import interval, inf, imath


On Turk, these packages installed in the usual way without fuss.

On bender, there was a problem with the how libcrlibm.a was build:
relocation R_X86_64_32S against `.rodata' can not be used when making
a shared object; recompile with -fPIC

To fix this, in the crlibm source tree,in scs_lib/Makefile add the -fPIC 
option in this line:
CPPFLAGS = -fPIC 


On a Mac, you need to do this:

I downloaded crlibm-1.0beta4.tar.gz from 
http://lipforge.ens-lyon.fr/www/crlibm/ 

$ untar
$ configure
$ make
$ sudo make install. 

Then I downloaded pyinterval-1.0b21.tar.gz from
http://pypi.python.org/pypi/pyinterval
and untarred it.


Here is the tricky part.

in setup.py pyinterval, change the lines near the end from this

             include_dirs = ['/opt/crlibm/include'],
             library_dirs = ['/opt/crlibm/lib'],

to this


            include_dirs = ['/usr/local/include'],
            library_dirs = ['/usr/local/lib'],

and then 

python setup.py build
sudo python setup.py install

The tests in cosmic/TestTimeMask.py uses the class Cosmic which uses
this, and it will also be used to locate times to mask hot pixels.


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

/hotpix contains tools for finding location and time of hot (and possibly 'cold') pixels

/photonlist contains tools for creating calibrated photon lists and creating stacked, RA/dec-mapped images from them.

Each directory contains a /test subdirectory, where code to test the module will be stored.

***

A beginner's guide and the observing log for the Palomar 2012 run are in examples/demos

***
This document uses the markdown syntax, see http://daringfireball.net/projects/markdown/


***
Getting it working on bender during the Palomar 2013 run.  

Bender has the canopy version of python installed.  To get pyqt working:

yum install qt4

I then tried yum install PyQt4 but the canopy version of python did
not pick it up.

sudo scp -r stoughto@turk.physics.ucsb.edu:/usr/local/epd/lib/python2.7/site-packages/PyQt4 .
This needs the same version of sip.  On turk, sip -V says it is version 4.13.2 so I downloaded sip-4.13.3.zip,
unzip, python configure.py; make; sudo make install.


