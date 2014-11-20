Hot pixel time-masking (and probably other types of masking soon too)
---------------------------------------------------------------------

1) To create a hot pixel time-mask file, use findHotPixels. E.g.:

		import hotPixels as hp
        hotpix.findHotPixels('obs_20121211-024511.h5', 'testOutput.h5', 
                      'hotPixels.dict', startTime=0, endTime=5)
        
        - Writes time-mask data for 'obs_2012...' to file 'testOutput.h5', 
            using parameter file 'hotPixels.dict', and considering only 
            the first 5 seconds of the input file (here the start/end time 
            values override those in the parameter file).


2) To use hot pixel time-masks, use ObsFile methods loadHotPixCalFile(),
	switchOnHotPixTimeMask(), switchOffHotPixTimeMask(), and getPixelBadTimes():
	
	.loadHotPixCalFile(hotpixfilename) - loads in a hot pixel time-mask file for
		the current ObsFile instance, and automatically switches it on. Hot pixel
		time masking should then automatically be applied to pretty much any
		relevant ObsFile calls (getPixelCountImage, getPixelSpectrum, etc.)
		without needing to do anything further.
	.switchOffHotPixTimeMask() - switches off hot pixel time masking (though
		leaves it still loaded into the current obsFile instance).
	.swichOnHotPixTimeMask() switches the masking back on. (Automatically called
		by 'loadHotPixelCalFile() unless otherwise requested).
	.getPixelBadTimes(pixelRow,pixelCol) - returns the time interval(s) for 
		which a given pixel is bad (hot or whatever else may come up in the
		future).
	
	
	E.g. to get a masked image:
	
		import ObsFile as of
		myObsFile = of.ObsFile('obs_20121209-044636.h5')
		myObsFile.loadHotPixCalFile('hotPix_20121209-044636.h5')
		image = myObsFile.getPixelCountImage(firstSec=3.5, integrationTime=9.2)
	
	- Gets an 9.2sec long image integration starting at second 3.5 into variable
	'image'.
	
	**NOTE - NON-INTEGER VALUES CAN NOW BE SUPPLIED FOR getTimedPacketList() AND
	THEREFORE FOR ALL CURRENT OBSFILE ROUTINES WHICH TAKE TAKING START AND
	INTEGRATION TIMES**
	
	

3) Other possibly helpful functions in hotPixels:
	
		readHotPixels: readHotPixels: reads in a hot-pixel .h5 file into somewhat 
				sensible structures for use by external routines.

		checkInterval: creates a 2D mask for a given time interval within a given 
                exposure.
                

4) Useful routines for diagnosing hot pixels:

    - hotpix/quantifyBadTime.py - spits out diagnostics on the occurence of hot 
                            pixels for a raw obs file.
    - In util/ObsFile.py:
        . plotPixelLightCurve() - (*NEW*) plots a light curve for a given pixel, marking
                            intervals where the pixel is flagged as bad (assuming
                            a hot pixel mask is already loaded).
    - Note that hotpix.findHotPixels() also has a bunch of options for displaying what
                            it's doing as it runs.



More documentation is included in the source files themselves.


- Updated, JvE 10/28/2014.
