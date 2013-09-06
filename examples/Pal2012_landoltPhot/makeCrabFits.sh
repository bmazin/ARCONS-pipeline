#! /bin/bash

for N in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
do

echo "Running Crab PSF fit $N "
python fitPsf.py crabparams.dict $N

done

./makeCrabSpectra.sh

echo "done"
