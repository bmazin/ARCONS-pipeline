import numpy as np
import string
def getPlotValues(hg, **kwargs):
    """
    return lists of x and y points for plotting a histogram in steps

    hg can be either the tuple generated by numpy.histogram or
    an iterable of bin population
    """
    ylog = kwargs.get('ylog',False)
    if isinstance(hg,tuple):
        x = hg[1]
        y = hg[0]
    else:
        x = range(len(hg)+1)
        y = hg


    yPlotMin = getYPlotMin(y,ylog)
    xp = []
    yp = [yPlotMin]
    for i in range(len(y)):
        yval = y[i]

        xval = x[i]
        xp.append(xval)
        xp.append(xval)
        yp.append(yval)
        yp.append(yval)

    yp.append(yPlotMin)
    xPlotMax = x[-1]
    xp.append(xPlotMax)
    xp.append(xPlotMax)
    return xp,yp

def getYPlotMin(y,ylog):
    """
    return the minimum of the y values, ignoring non-positive values 
    if ylog is true
    """
    retval = float('inf')
    for yVal in y:
        if not (ylog and yVal <= 0):
            retval = min(retval, yVal)
    return retval
