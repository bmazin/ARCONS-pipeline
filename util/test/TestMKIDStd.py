import unittest
import numpy as np
import math
from util import MKIDStd
class TestMKIDStd(unittest.TestCase):

    def testMeasureBandPassFlux(self):
        # taken from "Spectrophotometry and Magnitudes by Chris Churchill
        # table 5.1
        # http://ganymede.nmsu.edu/cwc/Teaching/ASTR605/Lectures/mags-extinct.pdf

        answers = { 'U':3.18e-9, 'B':6.60e-9, 'V':3.64e-9, \
                        'R':3.08e-9, 'I':2.55e-9}
        std = MKIDStd.MKIDStd()
        aFlux = std.load("vega")
        for filter in ['U','B','V','R','I']:
            aFilter = std.filters[filter]
            bpf = std.measureBandPassFlux(aFlux,aFilter)
            print "filter=",filter,"  answer=",answers[filter], "  bpf=",bpf

    def testGetVegaMag(self):
        """
        confirm that the magnitude of Vega in all bands is 0.03
        """
        std = MKIDStd.MKIDStd()
        vegaFlux = std.load("vega")
        bd17Flux = std.load("bd17")
        for filter in ['U','B','V','R','I']:
            aFilter = std.filters[filter]            
            mag = std.getVegaMag(vegaFlux, aFilter)
            self.assertAlmostEqual(0.03, mag, \
                                       msg="filter=%s mag=%f"%(filter,mag))

    def testCalspecMags(self):
        """
        compare values from getVegaMag for BD17 to values in
        http://www.stsci.edu/hst/observatory/cdbs/calspec.html
        """
        std = MKIDStd.MKIDStd()
        bFilter = std.filters['B']
        vFilter = std.filters['V']

        # BD17
        bd17Flux = std.load("bd17")
        B = std.getVegaMag(bd17Flux, bFilter)
        V = std.getVegaMag(bd17Flux, vFilter)
        self.assertAlmostEqual(B-V, 0.44, places=1, msg="value=%f"%B)
        self.assertAlmostEqual(B, 9.47, places=0, msg="value=%f"%B)

        
    def testConvert(self):
        """
        Test the equation for converting AB magnitude to flux and back again, 
        from http://en.wikipedia.org/wiki/AB_magnitude
        """
        for i in range(100):
            AB = 10 + i/10.0
            wave = 1000 + 100*i
            f = (10**(-2.406/2.5))*(10**(-0.4*AB))/(wave**2)
            AB_new = -2.5*math.log10(f) - 5*math.log10(wave) - 2.406
            #print "AB=",AB,"  wave=",wave, " f=",f, "  AB_new=",AB_new
            self.assertAlmostEqual\
                (AB, AB_new, msg="AB=%f f=%f AB_new=%f" % (AB, f, AB_new))


    def testGetIndex(self):
        """
        demonstrate how np.searchsorted works
        """
        w = 3000+2*np.arange(1000)
        referenceWavelength = 4567.89
        index = np.searchsorted(w, referenceWavelength);
        
        print "w[index-1]=",w[index-1]
        print "referenceWavelength=",referenceWavelength
        print "w[index]=",w[index]

        assertTrue(w[index-1] <= referenceWavelength)
        assertTrue(referenceWavelength <= w[index])

if __name__ == '__main__':
    unittest.main()



