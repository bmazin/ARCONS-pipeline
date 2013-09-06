#! /bin/bash

#for N in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 
for N in 30 31 32 33
do

echo "Running Crab file $N with "
python makeSpectralFrames.py fluxParams/fluxParams_crab_$N.txt $N

done

mv *.npz ~/scratch/standards/
mv *.gif ~/scratch/standards/

for N in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
do

echo "Running Crab PSF fit $N "
python fitPsf.py crabparams.dict $N

done


for N in 25 26 27 28 29 30 31 32 33
do

echo "Running Crab PSF fit $N "
python fitPsf.py crabNight2params.dict $N

done

./makeCrabSpectra.sh

echo "done"
