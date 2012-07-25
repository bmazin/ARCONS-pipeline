Description of what is in the examples folder of ARCONS-pipeline

git-cheat-sheet.txt -- a very terse explanation of the basics of git.
		    You need to have an account on github and Ben has
		    to give you access to ARCONS-pipeline for this to
		    work.

simple-xray -- contains a python file and one hdf5 file that reads all the
	    photon events and makes two histograms.


palomar-2011 -- demonstrates how to read files from the Palomar run in 2011.  

	     * palomar-2011-count-photons.py counts the total number of
	     photons in a file.  It uses /beammap/beamimage to
	     decide what to iterate.

	     * palomar-2011-to-fits.py Uses /beammap/beamimage to make a
	     2d histogram of the counts in each pixel, and writes as the
	     file new.fits
