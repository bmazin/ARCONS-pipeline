from flatcal.flatCal import FlatCal
from util.ObsFile import ObsFile
import matplotlib.pyplot as plt
import numpy as np

def main():
    np.set_printoptions(threshold=np.nan)
    #obs_20120919-131142.h5,obs_20120919-131346.h5
    ob = ObsFile('obs_20120919-131142.h5')
    ob.loadWvlCalFile('calsol_20120917-072537.h5')
    ob.loadFlatCalFile('flatsol_20120919-131142.h5')
    ob.createPhotonList('pl_20120919-131142.h5')


if __name__ == '__main__':
    main()
