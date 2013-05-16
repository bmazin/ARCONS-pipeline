import numpy as np

MJDtcbfile = np.loadtxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3MJD-TDB.txt')
MJDtdbfile = []
for line in MJDtcbfile:
    MJDtcb = line[1]
    MJDtdb = MJDtcb+(-1.550519768*10**-8*(MJDtcb+2400000.5-2443144.5003725)*86400-6.55*10**-5)/(24.*3600.)
    MJDtdbfile.append(MJDtdb)

np.savetxt('/Scratch/dataProcessing/SDSS_J0926/Dec8fitpsfBlueInt3longMJD_TDB.txt',MJDtdbfile,fmt='%.12f')
