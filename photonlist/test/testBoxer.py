import numpy as np
from photonlist.boxer import boxer

def testBoxer(i=0.0, j=0.0,
              x=np.array([-1.0, 0.0, 1.0, 0.0]),
              y=np.array([0.0, 1.0, 0.0, -1.0])):

    #Corners of a 2x2 square rotated at 45deg and centered at origin:
    rotSqX = np.array([-1.,0.,1.,0.])
    rotSqY = np.array([0.,1.,0.,-1.])
    #Integer values for centers of unit square:
    iList=[0,1,1,0,0,1,1,0]
    jList=[0,0,1,1,0,0,1,1]
    #List of quadrilaterals:
    xList=[rotSqX,rotSqX,rotSqX,rotSqX,
           rotSqX+0.5,rotSqX+0.5,rotSqX+0.5,rotSqX+0.5]
    yList=[rotSqY,rotSqY,rotSqY,rotSqY,
           rotSqY,rotSqY,rotSqY,rotSqY]
    #List of expected overlap areas:
    expectedAreaList=[1.0,0.25,0,0.25,0.75,0.75,0.125,0.125]

    print 'i,j,x,y,calculatedArea,expectedArea'
    for i,j,x,y,expectedArea in zip(iList,jList,xList,yList,expectedAreaList):
        calculatedArea=boxer(i, j, x, y)
        print 'i,j: ',i,j
        print 'x: ',x
        print 'y: ',y
        print 'Calculated area: ',calculatedArea
        print 'Expected area: ',expectedArea
        print
        assert np.abs(calculatedArea-expectedArea) < 1e-7
        

    print 'All good.'