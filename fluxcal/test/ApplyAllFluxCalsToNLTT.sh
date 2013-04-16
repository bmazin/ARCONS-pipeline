#! /bin/bash

for FILE1 in fluxParams_nltt11748.txt 
do

for FILE2 in fluxParams_hiltner600.txt fluxParams_bd25.txt fluxParams_hd849a.txt
do
echo "Running Cal on $FILE1 with $FILE2"
python testApertureSpectrum.py $FILE1 $FILE2 0
done

for FILE2 in fluxParams_bd25_1.txt fluxParams_hd849a_1.txt
do
echo "Running Cal on $FILE1 with $FILE2"
python testApertureSpectrum.py $FILE1 $FILE2 1
done

done

