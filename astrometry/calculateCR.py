#!/bin/python

'''
Author: Paul Szypryt		Date: April 30, 2013

Based on CircleFit.java (September 22, 2011 Jennifer Milburn, Original Implementation).

Takes a list of x and y positions on the array and calculates the location of a circle center that
would fit those points.
'''

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm

times, x, y, ha =np.loadtxt('/home/pszypryt/Scratch/centroid_test/centroid_list.txt',unpack='True',skiprows=1)  
numberOfPoints = len(x)


x=np.array([37.7630968422, 37.59877992885, 37.236499667099999, 37.177924051299996, 36.783136532100002, 36.271422505549999, 35.78070321165, 35.442410522149999, 35.037066280399998, 34.461779605399997, 33.566525930650002, 32.026794121449996, 31.579013769550002])
y=np.array([19.826419748699998, 20.511769704350002, 21.3182756912, 21.878565210700003, 22.35456673885, 23.730502243049997, 24.615149849150001, 24.944072180550002, 25.500067872800003, 26.695203749000001, 27.6978730437, 29.844366393850002, 30.4193377876])

meanX = x.mean()
meanY = y.mean()

xValues = x-meanX
yValues = y-meanY

sU = xValues
sUU = xValues * xValues
sUUU = xValues * xValues * xValues
sV = yValues
sVV = yValues * yValues
sVVV = yValues * yValues * yValues
sUV = xValues*yValues
sUVV = xValues * yValues * yValues
sVUU = yValues * xValues * xValues

sumU = np.sum(sU)
sumUU = np.sum(sUU)
sumUUU = np.sum(sUUU)
sumV = np.sum(sV)
sumVV = np.sum(sVV)
sumVVV = np.sum(sVVV)
sumUV = np.sum(sUV)
sumUVV = np.sum(sUVV)
sumVUU = np.sum(sVUU)

print 'Su = ' + str(sumU)
print 'Suu = ' + str(sumUU)
print 'Suuu = ' + str(sumUUU)
print 'Sv = ' + str(sumV)
print 'Svv = ' + str(sumVV)
print 'Svvv = ' + str(sumVVV)
print 'Suv = ' + str(sumUV)
print 'Suvv = ' + str(sumUVV)
print 'Svuu = ' + str(sumVUU)
print 'MeanX = ' + str(meanX)
print 'MeanY = ' + str(meanY)

centerX = -(((-sumUV*sumVUU)+(sumUUU*sumVV)+(sumUVV*sumVV)+(-sumUV*sumVVV))/(2.0*((sumUV*sumUV) -sumUU*sumVV)))
centerY = -(((-sumUUU*sumUV)+(-sumUV*sumUVV)+(sumUU*sumVUU)+(sumUU*sumVVV))/(2.0*((sumUV*sumUV) -sumUU*sumVV)))

#print 'CenterX = ' + str(centerX)
#print 'CenterY = ' + str(centerY)

alpha = (centerX*centerX)+(centerY*centerY)+((sumUU+sumVV)/numberOfPoints)

print 'Alpha = ' + str(alpha)

centerX = centerX + meanX
centerY = centerY + meanY
radius = np.sqrt(alpha)
print 'CenterX = ' + str(centerX) + ' CenterY = ' + str(centerY) + ' Radius = ' + str(radius)

fig = plt.figure()
circle=plt.Circle((centerX,centerY),radius,fill=False)
fig.gca().add_artist(circle)
'''
xFitP = np.linspace(centerX-radius,centerX+radius,500)
xFitN = xFitP[::-1]
yFitP = centerY + np.sqrt(radius*radius-(xFitP-centerX)*(xFitP-centerX))
yFitN = centerY - np.sqrt(radius*radius-(xFitN-centerX)*(xFitN-centerX))
xFit = np.append(xFitP,xFitN)
yFit = np.append(yFitP,yFitN)
'''
plt.plot(x,y,'b.', centerX,centerY,'ro')
plt.xlim([0,43])
plt.ylim([0,45])
plt.show()
