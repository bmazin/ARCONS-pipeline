import numpy as np

#This is used to convert the .txt file of barycentered MJD(TCB) from Tempo 2 to MJD(TDB).

MJDtcbfile = np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAll.MJDlong-TCB.txt')
MJDtdbfile = []
for line in MJDtcbfile:
    MJDtcb = line[1]
    MJDtdb = MJDtcb+(-1.550519768*10**-8*(MJDtcb+2400000.5-2443144.5003725)*86400-6.55*10**-5)/(24.*3600.)
    MJDtdbfile.append(MJDtdb)

np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfAllInt.1MJD_TDB.txt',MJDtdbfile,fmt='%.12f')
