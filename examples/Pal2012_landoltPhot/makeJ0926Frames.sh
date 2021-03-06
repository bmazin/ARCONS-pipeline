#! /bin/bash

for N in 0 1 2 3 4 5 6 7 8 9 10 11 12
do

echo "Running J0926 file $N"
python makeSpectralFrames.py fluxParams/fluxParams_sdssj0926_$N.txt $N

mv *.npz ~/scratch/standards
mv *.gif ~/scratch/standards

python fitPsf.py sdssj0926params.dict $N

done
echo "done"
