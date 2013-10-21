'''
Author: Julian van Eyken                Date: May 6 2013

Any global constants for general use in the pipeline. Only things that
can't be obtained from file headers, and things that are genuinely constant
and not likely to change, or at least VERY rarely. Things that are not
suitable as user input parameters, and nothing that is specific to a 
particular instrument or instrument configuration. Any weird or unusual
physical constants (that aren't already in the astropy.constants package)
can go here too.
'''

#Define as globals. E.g.:

v_thresh = 88    #Velocity threshold for time travel (for DeLorean, mph)
P_crit = 1.21    #Flux capacitor critical power level (GW)
