from matplotlib import rcParams, rc
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

# common setup for matplotlib
params = {'savefig.dpi': 300, # save figures to 300 dpi
          'axes.labelsize': 14,
          'lines.linewidth': 1.5,
          'text.fontsize': 14,
          'legend.fontsize': 14,
          'xtick.labelsize': 14,
          'ytick.major.pad': 6,
          'xtick.major.pad': 6,
          'ytick.labelsize': 14}
# use of Sans Serif also in math mode
rc('text.latex', preamble='\usepackage{sfmath}')

rcParams.update(params)

p1 = np.load('sigP1.npz')
p2 = np.load('sigP2.npz')

idxOffsets=p1['idxOffsets']
p1Sig = p1['nSigmaByIdxOffset']
p2Sig = p2['nSigmaByIdxOffset']


fig = plt.figure(figsize=(8.,6.))
ax = fig.add_axes([.1,.4,.85,.5])
ax2 = fig.add_axes([.1,.1,.85,.3])

ax.plot(idxOffsets,np.abs(p2Sig),'k')
fig.text(.04,.9,'Excess Peak Optical Emission (Standard Deviations)',size=14,rotation='vertical')
#ax2.set_ylabel('Standard Deviations of Peak Height from Average Peak')
ax2.set_xlabel('Pulse Offset Relative to GRP (number of periods)')
ax.set_ylim((0,8))
ax2.set_ylim((0,5))


ax2.plot(idxOffsets,np.abs(p1Sig),'k')

ax.xaxis.set_visible(False)
ax.xaxis.set_ticks([])

fig.text(.125,.85,'(a)',size=16)
ax2.yaxis.get_major_ticks()[-1].label1.set_visible(False)
fig.text(.125,.35,'(b)',size=16)

plt.show()
