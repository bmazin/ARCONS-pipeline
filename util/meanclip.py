import numpy
from util import robust_sigma

def meanclip(indata, clipsig=3.0, maxiter=5, converge_num=0.02, verbose=0):
   """
   Computes an iteratively sigma-clipped mean on a
   data set. Clipping is done about median, but mean
   is returned.

   .. note:: MYMEANCLIP routine from ACS library.

   :History:
       * 21/10/1998 Written by RSH, RITSS
       * 20/01/1999 Added SUBS, fixed misplaced paren on float call, improved doc. RSH
       * 24/11/2009 Converted to Python. PLL.
       * 18/10/2012 shamelessly stolen from https://gist.github.com/1310947 by Chris S.
       * 04/12/2012 Modified to import robust_sigma and return nSurvived
   Examples
   --------
   >>> mean, sigma, nSurvived = meanclip(indata)

   Parameters
   ----------
   indata: array_like
       Input data.

   clipsig: float
       Number of sigma at which to clip.

   maxiter: int
       Ceiling on number of clipping iterations.

   converge_num: float
       If the proportion of rejected pixels is less than
       this fraction, the iterations stop.

   verbose: {0, 1}
       Print messages to screen?

   Returns
   -------
   mean: float
       N-sigma clipped mean.

   sigma: float
       Standard deviation of remaining pixels.

   """
   # Flatten array
   skpix = indata.reshape( indata.size, )

   ct = indata.size
   iter = 0; c1 = 1.0 ; c2 = 0.0

   while (c1 >= c2) and (iter < maxiter):
       lastct = ct
       medval = numpy.median(skpix)
       sig = numpy.std(skpix)
       wsm = numpy.where( abs(skpix-medval) < clipsig*sig )
       ct = len(wsm[0])
       if ct > 0:
           skpix = skpix[wsm]

       c1 = abs(ct - lastct)
       c2 = converge_num * lastct
       iter += 1
   # End of while loop

   mean  = numpy.mean( skpix )
   sigma = robust_sigma.robust_sigma( skpix )

   if verbose:
       prf = 'MEANCLIP:'
       print '%s %.1f-sigma clipped mean' % (prf, clipsig)
       print '%s Mean computed in %i iterations' % (prf, iter)
       print '%s Mean = %.6f, sigma = %.6f' % (prf, mean, sigma)

   nSurvived = len(wsm[0])
   return mean, sigma, nSurvived


