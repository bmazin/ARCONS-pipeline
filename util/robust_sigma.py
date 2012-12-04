import numpy

def robust_sigma(in_y, zero=0):
   """
   Calculate a resistant estimate of the dispersion of
   a distribution. For an uncontaminated distribution,
   this is identical to the standard deviation.

   Use the median absolute deviation as the initial
   estimate, then weight points using Tukey Biweight.
   See, for example, Understanding Robust and
   Exploratory Data Analysis, by Hoaglin, Mosteller
   and Tukey, John Wiley and Sons, 1983.

   .. note:: ROBUST_SIGMA routine from IDL ASTROLIB.

   :History:
       * H Freudenreich, STX, 8/90
       * Replace MED call with MEDIAN(/EVEN), W. Landsman, December 2001
       * Converted to Python by P. L. Lim, 11/2009

   Examples
   --------
   >>> result = robust_sigma(in_y, zero=1)

   Parameters
   ----------
   in_y: array_like
       Vector of quantity for which the dispersion is
       to be calculated

   zero: int
       If set, the dispersion is calculated w.r.t. 0.0
       rather than the central value of the vector. If
       Y is a vector of residuals, this should be set.

   Returns
   -------
   out_val: float
       Dispersion value. If failed, returns -1.

   """
   # Flatten array
   y = in_y.reshape(in_y.size, )

   eps = 1.0E-20
   c1 = 0.6745
   c2 = 0.80
   c3 = 6.0
   c4 = 5.0
   c_err = -1.0
   min_points = 3

   if zero:
       y0 = 0.0
   else:
       y0 = numpy.median(y)

   dy    = y - y0
   del_y = abs( dy )

   # First, the median absolute deviation MAD about the median:

   mad = numpy.median( del_y ) / c1

   # If the MAD=0, try the MEAN absolute deviation:
   if mad < eps:
       mad = numpy.mean( del_y ) / c2
   if mad < eps:
       return 0.0

   # Now the biweighted value:
   u  = dy / (c3 * mad)
   uu = u*u
   q  = numpy.where(uu <= 1.0)
   count = len(q[0])
   if count < min_points:
       print 'ROBUST_SIGMA: This distribution is TOO WEIRD! Returning', c_err
       return c_err

   numerator = numpy.sum( (y[q]-y0)**2.0 * (1.0-uu[q])**4.0 )
   n    = y.size
   den1 = numpy.sum( (1.0-uu[q]) * (1.0-c4*uu[q]) )
   siggma = n * numerator / ( den1 * (den1 - 1.0) )

   if siggma > 0:
       out_val = numpy.sqrt( siggma )
   else:
       out_val = 0.0

   return out_val
