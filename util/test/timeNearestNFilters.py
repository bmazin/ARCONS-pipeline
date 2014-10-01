'''
To set up a timing test on the 'nearest-N' filters in util.utils, since it
runs pretty slowly at the moment.
'''

import os.path
from util import utils
import timeit
import numpy as np

testImage = os.path.join(os.path.dirname(__file__), 'egImageArray.npy')


print 'Timing: utils.nearestNrobustSigmaFilter()'
t1 = timeit.timeit('utils.nearestNrobustSigmaFilter(im)',
        setup='from util import utils; import numpy; im = numpy.load("' + testImage + '")',
        number=3)
print 'Done - average '+str(t1/3.)+'sec per call'

print 'Timing: utils.nearestNmedFilter()'
t2 = timeit.timeit('utils.nearestNmedFilter(im)',
        setup='from util import utils; import numpy; im = numpy.load("' + testImage + '")',
        number=3)
print 'Done - average '+str(t2/3.)+'sec per call'

print 'Timing: utils.findNearestFinite()'
t3 = timeit.timeit('utils.findNearestFinite(im,10,10,n=25)',
        setup='from util import utils; import numpy; im = numpy.load("' + testImage + '")',
        number=1000)
print 'Done - average '+str(t3/1000.)+'sec per call'

