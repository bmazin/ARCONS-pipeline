import unittest
import numpy as np
import math
from util import MKIDStd
class TestMKIDStd(unittest.TestCase):

    def testPlot(self):
        std = MKIDStd.MKIDStd()
        std.plot()

    def testGetVegaMag(self):
        """
        confirm that the magnitude of Vega in all bands is 0.03
        """
        std = MKIDStd.MKIDStd()
        vegaFlux = std.load("vega")
        for filter in ['U','B','V','R','I']:
            mag = std.getVegaMag("vega", filter)
            self.assertAlmostEqual(0.03, mag, \
                                       msg="filter=%s mag=%f"%(filter,mag))

    def testCalspecMags(self):
        """
        compare values from getVegaMag for BD17 to values in
        http://www.stsci.edu/hst/observatory/cdbs/calspec.html
        """
        std = MKIDStd.MKIDStd()

        # BD17
        B = std.getVegaMag("bd17", "B")
        V = std.getVegaMag("bd17", "V")
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
        

        self.assertTrue(w[index-1] <= referenceWavelength)
        self.assertTrue(referenceWavelength <= w[index])

if __name__ == '__main__':
    unittest.main()



