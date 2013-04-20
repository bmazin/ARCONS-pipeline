#!/bin/python

'''
Author: Paul Szypryt		Date: March 8, 2013

Uses raw photon data (no fits, too small timescales) to calculate the Lomb-Scargle Periodogram.  Does this in a specified
aperture of given x and y location and radius.  Applies wavelength cals, flat cals, and other typical ObsFiles functions.
Averages over selected timesteps, optimized for high frequencies (>1 Hz).
'''
