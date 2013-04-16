#! /bin/bash

for FILE1 in fluxParams_hiltner600.txt fluxParams_HR3454.txt fluxParams_HR9087.txt fluxParams_feige66.txt fluxParams_bd25.txt fluxParams_hd849a.txt fluxParams_G158-100.txt
do

for FILE2 in fluxParams_hiltner600.txt fluxParams_HR3454.txt fluxParams_HR9087.txt fluxParams_feige66.txt fluxParams_bd25.txt fluxParams_hd849a.txt
do

echo "Running Cal on $FILE1 with $FILE2"
python testFluxCalApplication.py $FILE1 $FILE2

done
done
