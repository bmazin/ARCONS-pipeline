# A simple demonstration of how to use ephem to deal with segesimal strings

import ephem
import math

raSexigesimal = "15:30:00"
decSexigesimal = "-10:01:00"

ra = ephem.hours(raSexigesimal) 
dec = ephem.degrees(decSexigesimal) # oddly enough, "degrees" is in radians

print "ra=",ra
print "dec=",dec

print "ra in decimal degrees is ",ra.real*12/math.pi
print "dec in decimal degrees is ",dec.real*180/math.pi

