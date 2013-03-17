#! /bin/bash

for FILE1 in fluxParams_G158-100.txt fluxParams_crab.txt fluxParams_crabNight2.txt fluxParams_nltt11748.txt fluxParams_nltt11748Night2.txt fluxParams_sdssj0926.txt 
do

for FILE2 in fluxParams_hiltner600.txt
do
echo "Running Cal on $FILE1 with $FILE2"
python testApertureSpectrumNoFlux.py $FILE1 $FILE2 0
done

done

for FILE1 in fluxParams_G158-100_1.txt fluxParams_crab_1.txt fluxParams_crabNight2_1.txt fluxParams_sdssj0926_1.txt
do

for FILE2 in fluxParams_hiltner600.txt
do
echo "Running Cal on $FILE1 with $FILE2"
python testApertureSpectrumNoFlux.py $FILE1 $FILE2 2
done

done

