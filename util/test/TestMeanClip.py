import unittest
from util import meanclip
import numpy as np

class TestMeanclip(unittest.TestCase):
    """
    Generate 1000 gaussian normal numbers with three outliers.  Check
    that meanclip does the right thing
    """
    def testSimple(self):
        mu = 100
        sig = 10
        s = np.random.normal(mu, sig, 1000)
        s[400] = mu+13*sig
        s[500] = mu+15*sig
        s[654] = mu-15*sig
        mean,sigma,nSurvived = meanclip.meanclip(s, clipsig=9.0)
        self.assertEqual(nSurvived, 997, \
                             "three should get cut nSurived=%d" % nSurvived)
        self.assertAlmostEqual(mu, mean, delta=0.4, \
                                   msg="mu=%f mean=%f" % (mu,mean))
        self.assertAlmostEqual(sig, sigma, delta=0.6, \
                                   msg="sig=%f sigma=%f" % (sig,sigma))
if __name__ == '__main__':
    unittest.main()
