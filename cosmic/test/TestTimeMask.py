import os
import unittest
from headers import TimeMask
import tables
from tables.nodes import filenode
from util import FileName
from cosmic.Cosmic import Cosmic
import numpy as np
import matplotlib.pyplot as plt
import inspect
class TestTimeMask(unittest.TestCase):
    """
    test the TimeMask class and practice writing hdf5 files
    """


    def testEnum(self):
        """
        Write four rows into file enum.h5 and read it back.
        Note that the actual enum is stored inside the file
        """
        h5f = tables.openFile('enum.h5', 'w')

        fnode = filenode.newNode(h5f, where='/', name='timeMaskHdr')
        fnode.attrs.beginTime = 1234
        fnode.attrs.endTime = 5678
        fnode.close()

        tbl =  h5f.createTable('/','timeMask',TimeMask.TimeMask,"Time Mask")

        b = [1010000, 1020000, 1010001, 1030000]
        e = [1020000, 1030000, 1020001, 1040000]
        r = ["Flash in r0", "Flash in r2", "Flash in r7", "Merged Flash"]
        for i in range(len(b)):
            row = tbl.row
            row['tBegin'] = b[i]
            row['tEnd']   = e[i]
            row['reason'] = TimeMask.timeMaskReason[r[i]]
            row.append()
            tbl.flush()
        tbl.close()
        h5f.close()

        fid = tables.openFile("enum.h5", mode='r')

        node = fid.root.timeMaskHdr
        self.assertEqual(node.attrs.beginTime, 1234)
        self.assertEqual(node.attrs.endTime, 5678)

        table = fid.getNode('/timeMask')
        self.assertEqual(len(b), table.nrows)
        enum = table.getEnum('reason')
        for i in range(table.nrows):
            self.assertEqual(table[i]['tBegin'], b[i])
            self.assertEqual(table[i]['tEnd'], e[i])
            self.assertEqual(enum(table[i]['reason']),r[i])
        fid.close()

    def testCosmic(self):
        run = 'PAL2012'
        sundownDate = '20121211'
        obsDate = '20121212'
        seq = '121229'
        fn = FileName.FileName(run, sundownDate, obsDate+"-"+seq)
        cosmic = Cosmic(fn, beginTime=123, endTime=125)

        dictionary0 = cosmic.findCosmics(threshold=500, populationMax=2000)
        print "dictionary0 keyst=",dictionary0.keys()
        print "interval=",dictionary0['interval']

        hist = dictionary0['populationHg'][0]
        bins = np.arange(len(hist))
        plt.plot(bins, hist, label="no mask")
        
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.savefig(inspect.stack()[0][3]+".png")

if __name__ == '__main__':
    unittest.main()
