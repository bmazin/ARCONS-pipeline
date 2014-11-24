from mpltools import style
style.use('jfm')

import matplotlib as mpl

mpl.rcParams['figure.figsize'] = [5.0,3.3]
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['savefig.dpi'] = 300
#mpl.rcParams['savefig.bbox'] = 'tight'
#mpl.rcParams['savefig.pad_inches'] = 0.05
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = ['Computer Modern Roman']
mpl.rcParams['text.usetex'] = True
