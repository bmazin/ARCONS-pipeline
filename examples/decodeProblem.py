import tables

fullFileName = "/ScienceData/LICK2012/20120919/obs_20120920-092626.h5"
pixel = "/r2/p95/t1348133188"

fid = tables.openFile(fullFileName, mode='r')

secs = fid.getNode(pixel)

sec = secs[0]

timeMask = int(20*'1',2)   #bitmask of 20 ones

time = sec & timeMask

print "time = %d" % time[0]

sec0 = sec[0]

sec

sec0

print "this will not work"
masked = sec0 & timeMask


