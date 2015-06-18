import os
import math
import unittest
from util import ObsFileSeq, FileName, ObsFile
import matplotlib.pyplot as plt
import matplotlib.patches as patches
class TestObsFileSeq(unittest.TestCase):

    def testGetScaleTheta(self):
        flip = -1
        scale = 0.4/3600 # degrees/pixel
        theta = 31.2 # degrees
        st = math.sin(math.radians(theta))
        ct = math.cos(math.radians(theta))

        rMax, cMax = 300, 300
        rStar, cStar = 140, 173
        fig, axes = plt.subplots()
        # draw row,col arrows
        axes.arrow(0, 0, cMax, 0, head_width=5, head_length=10, fc='k', ec='k')
        axes.arrow(0, 0, 0, rMax, head_width=5, head_length=10, fc='k', ec='k')
        rdl = 300

        # draw the ra, dec axes
        axes.arrow(cStar-flip*rdl*ct/2, rStar-flip*rdl*st/2, flip*rdl*ct, flip*rdl*st, head_width=5, head_length=10)
        axes.text(cStar-flip*rdl*ct/2+flip*rdl*ct, rStar-flip*rdl*st/2+flip*rdl*st-5, "ra",
                  size='smaller', va='top', ha='center')
        axes.arrow(cStar+rdl*st/2, rStar-rdl*ct/2, -rdl*st, rdl*ct, head_width=5, head_length=10)
        axes.text(cStar+rdl*st/2-rdl*st, rStar-rdl*ct/2+rdl*ct+5, "dec",
                  size='smaller', va='bottom', ha='center')
        axes.plot(cStar, rStar, 'r*')


        # define how large a frame is
        nRow, nCol = 100, 120
        # define where the star shows up on frame a
        rowa, cola = 80, 60
        # calculate r,c position of row,col = 0,0 in frame a
        row0a = rStar - rowa
        col0a = cStar - cola
        # draw frame a in blue
        axes.add_patch(patches.Rectangle((col0a, row0a), nCol, nRow, alpha=0.1, facecolor='blue'))

        # define where the star shows up on frame b
        rowb, colb = 50, 40
        # calculate r,c position of row,col = 0,0 in frameb
        row0b = rStar - rowb
        col0b = cStar - colb
        # draw frame b in blue
        axes.add_patch(patches.Rectangle((col0b, row0b), nCol, nRow, alpha=0.1, facecolor='red'))

        # calculate ra, dec for the (0,0) corner of frame a and b
        ra0a = (scale/flip) * (-cola*ct - rowa*st)
        dec0a = scale * (cola*st - rowa*ct)
        ra0b = (scale/flip) * (-colb*ct - rowb*st)
        dec0b = scale * (colb*st - rowb*ct)

        # draw the ra0a line
        cBeg = col0a
        rBeg = row0a
        raEnd = ra0a
        decEnd = 0
        cEnd = (flip*raEnd*ct - decEnd*st)/scale + col0a + cola
        rEnd = (flip*raEnd*st + decEnd*ct)/scale + row0a + rowa
        axes.plot((cBeg,cEnd), (rBeg, rEnd), "g:")
        axes.text(cEnd, rEnd, "ra0a", size='smaller', rotation=theta-90, ha='right', va='bottom', color='b')
        # draw the dec0a line
        cBeg = col0a
        rBeg = row0a
        raEnd = 0
        decEnd = dec0a
        cEnd = (flip*raEnd*ct - decEnd*st)/scale + col0a + cola
        rEnd = (flip*raEnd*st + decEnd*ct)/scale + row0a + rowa
        axes.plot((cBeg,cEnd), (rBeg, rEnd), "b:")
        axes.text(cEnd, rEnd, "dec0a", size='smaller', rotation=theta, ha='left', va='bottom', color='b')

        # draw the ra0b line
        cBeg = col0b
        rBeg = row0b
        raEnd = ra0b
        decEnd = 0
        cEnd = (flip*raEnd*ct - decEnd*st)/scale + col0b + colb
        rEnd = (flip*raEnd*st + decEnd*ct)/scale + row0b + rowb
        axes.plot((cBeg,cEnd), (rBeg, rEnd), "r:")
        axes.text(cEnd, rEnd, "ra0b", size='smaller', rotation=theta-90, ha='right', va='bottom', color='r')

        # draw the dec0b line
        cBeg = col0b
        rBeg = row0b
        raEnd = 0
        decEnd = dec0b
        cEnd = (flip*raEnd*ct - decEnd*st)/scale + col0b + colb
        rEnd = (flip*raEnd*st + decEnd*ct)/scale + row0b + rowb
        axes.plot((cBeg,cEnd), (rBeg, rEnd), "r:")

        # draw the rowa, cola line
        axes.plot((0, col0a), (row0a, row0a), "b--")
        axes.text(col0a, -2, "c0a", size='smaller', rotation=90, ha='center', va='top', color='b')
        axes.plot((col0a, col0a), (0, row0a), "b--")
        axes.text(-2, row0a, "r0a", size='smaller', ha='right', va='center', color='b')
        axes.plot((0, col0b), (row0b, row0b), "r--")
        axes.text(col0b, -2, "c0b", size='smaller', rotation=90, ha='center', va='top', color='r')
        axes.plot((col0b, col0b), (0, row0b), "r--")
        axes.text(-2, row0b, "r0b", size='smaller', ha='right', va='center', color='r')

        axes.set_xlim((-30, cMax+10))
        axes.set_ylim((-30, rMax+10))
        axes.set_aspect("equal")

        axes.set_xlabel('pixel column')
        axes.set_ylabel('pixel row')
        axes.set_title("flip=%d"%flip)
        plt.savefig("testGetScaleTheta.png")

        matchList = [
            dict(ra=ra0a, dec=dec0a, row=row0a, col=col0a),
            dict(ra=ra0b, dec=dec0b, row=row0b, col=col0b),
        ]
        scaleTheta = ObsFileSeq.ObsFileSeq.getScaleTheta(
                        matchList, flip=flip)
        print "scaleTheta=",scaleTheta
        print "scale=", scaleTheta['scale']*3600, " arcsec/pixel"
        print "theta=", math.degrees(scaleTheta['theta'])," degrees"

    def aatestGetScaleThetaFromMatches01(self):
        for flip in (1, -1):
            for theta in [-10.0, 30.0, -123.4, 98.7, 0.0, 90.0]:
                ct = math.cos(math.radians(theta))
                st = math.sin(math.radians(theta))
                for scale in (0.4/3600, 0.321/3600):
                    msg = "flip=%d theta=%f scale=%f" % \
                        (flip, theta, scale*3600)
                    dra, ddec = 1, 20  # in arcseconds
                    ra0a = 0
                    dec0a = 0
                    ra0b = dra/3600.0
                    dec0b = ddec/3600.0
                    cola = (flip*ra0a*ct - dec0a*st)/scale
                    rowa = (flip*ra0a*st + dec0a*ct)/scale
                    colb = (flip*ra0b*ct - dec0b*st)/scale
                    rowb = (flip*ra0b*st + dec0b*ct)/scale
                    matchList = [
                        {"ra": ra0a, "dec": dec0a, "row": rowa, "col": cola},
                        {"ra": ra0b, "dec": dec0b, "row": rowb, "col": colb},
                        ]
                    scaleTheta = ObsFileSeq.ObsFileSeq.getScaleTheta(
                        matchList, flip=flip)
                    self.assertAlmostEqual(theta,
                                           math.degrees(scaleTheta['theta']),
                                            msg=msg)
                    self.assertAlmostEqual(scale, scaleTheta['scale'],
                                           delta=1e-3/3600, msg=msg)
     
       
    def aaatestGetScaleThetaFromMatches(self):
        scale = 0.4/3600.0  # degrees/pixel
        theta = 20 # degrees
        ct = math.cos(math.radians(theta))
        st = math.sin(math.radians(theta))
        
        dra,ddec = 1,20
        
        ra0a = 0
        dec0a = 0
        ra0b =  dra/3600.0
        dec0b = ddec/3600.0
        cola = ( ra0a*ct - dec0a*st)/scale
        rowa = ( ra0a*st + dec0a*ct)/scale
        colb = ( ra0b*ct - dec0b*st)/scale
        rowb = ( ra0b*st + dec0b*ct)/scale
        matchList = [
            {"ra":ra0a, "dec":dec0a, "row":rowa, "col":cola},
            {"ra":ra0b, "dec":dec0b, "row":rowb, "col":colb},
            ]
        st = ObsFileSeq.ObsFileSeq.getScaleTheta(matchList)
        print math.degrees(st['theta'])
        print 3600.0*st['scale']
    def aatestObsFileSeq(self):
        name = 'ring-20141020'
        run = "PAL2014"
        date = "20141020"
        tsl = [
            #'20141021-033954',
            #'20141021-034532',
            '20141021-035035',
            '20141021-035538',
            #'20141021-040041',
            #'20141021-040544',
            #'20141021-041047',
            ]
        dt = 200
        ofs = ObsFileSeq.ObsFileSeq(name,run,date,tsl,dt)
        ofs.getFrameList()
        ofs.loadSpectralCubes()
        
        #sc = ofs.getSpectralCubeByFrame(0)
        #sc = ofs.getSpectralCubeByFrame(1)
if __name__ == '__main__':
    unittest.main()
