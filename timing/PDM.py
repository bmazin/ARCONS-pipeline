import numpy as np
from scipy import stats

def PDM(times, fluxes, frequencies, numberOfBins = 10, binWidth = 0.1):
    """
    Perform phase dispersion minimization.  Need to add option for flux errors.
    """

    # Offset time array to make t0 = 0
    zeroPoint = times[0]
    times -= zeroPoint

    # Total number of data points, frequencies
    numberOfData = len(times)
    numberOfFrequencies = len(frequencies)

    # Calculate width used to center the bins
    widthForCenter = 1/float(numberOfBins)

    dispersions = np.zeros(len(frequencies))

    # Loop through total number of frequencies
    for iFrequency in range(numberOfFrequencies):

        # Initialize array for number of points in bin, may need to place in loop
        numPoints = np.zeros(numberOfBins)
        binVariance = np.zeros(numberOfBins)

        # Convert times to phase folded on frequencies[iFrequency], sort times, fluxes, fluxErrors
        sortedPhases, sortedFluxes = FoldTimes(times, fluxes, frequencies[iFrequency])

        overallVariance = stats.tvar(sortedFluxes)

        # Loop through total number of bins
        for iBin in range(numberOfBins):         

            # Use 'binWidth' to determine the min/max values of the bin
            binCenter = (iBin+1)*widthForCenter - 0.5*widthForCenter
            binMin = binCenter - 0.5*binWidth
            binMax = binCenter + 0.5*binWidth

            # Pick out fluxes that have associated phase between binMin and binMax
            # Account for bins with phases < 0 and > 1
            sample = sortedFluxes[np.where(np.logical_or(np.logical_or(np.logical_and(sortedPhases < binMax, sortedPhases >= binMin),np.logical_and(sortedPhases - 1 < binMax, sortedPhases - 1 >= binMin)),np.logical_and(sortedPhases + 1 < binMax, sortedPhases + 1 >= binMin)))]
            numPoints[iBin] = len(sample)

            # Calculate the variances of individual bins
            if numPoints[iBin] > 1:
                binVariance[iBin] = stats.tvar(sample)
            else:
                binVariance[iBin] = 0.
 
        # Calculate overall variance for samples
        numerator = 0.
        denominator = 0.    
        for iBin in range(numberOfBins):
            numerator += (float(numPoints[iBin])-1)*binVariance[iBin]
            denominator += float(numPoints[iBin])
        denominator -= numberOfBins        
        sampleVariance = numerator/denominator

        # Calculate dispersion measure
        dispersions[iFrequency] = sampleVariance/overallVariance

    return dispersions
        

def FoldTimes(times, yData, frequency):
    """
    Subroutine used by PDM to fold the lightcurve by the given trial periods
    """
    # Convert time to phase
    period = 1/float(frequency)
    phases = times/period - (times/period).astype('int')

    # Sort arrays in ascending phase order
    sortedIndices = phases.argsort()
    sortedPhases = phases[sortedIndices]
    sortedYData = yData[sortedIndices]
    #sortedYErrors = yErrors[sortedIndices]
    
    return sortedPhases, sortedYData
    #return sortedPhases, sortedYData, sortedYErrors

