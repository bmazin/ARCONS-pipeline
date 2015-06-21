"""
Make all the FITS files and the plot location file.
Use this to find the information needed fo ofs2.py
"""
from util.ObsFileSeq import ObsFileSeq
from util.readDict import readDict

d = readDict()
d.read_from_file("mosaic.dict")
name = 'PSN234416a'
dt = 200
ofs = ObsFileSeq(name, d['run'], d['date'], d['timeStampList'], dt)
#ofs.plotLocations(name+".png")
wvMin = 3000
wvMax = 13000
ofs.makeAllFitsFiles(wvMin, wvMax)


