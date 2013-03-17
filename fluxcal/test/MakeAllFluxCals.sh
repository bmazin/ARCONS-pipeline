#! /bin/bash

for FILE1 in fluxParams_bd25_1.txt fluxParams_bd25.txt fluxParams_feige66_2.txt fluxParams_feige66.txt fluxParams_G158-100_1.txt fluxParams_G158-100_2.txt fluxParams_G158-100.txt fluxParams_hd849a_1.txt fluxParams_hd849a.txt fluxParams_hiltner600.txt fluxParams_HR3454_1.txt fluxParams_HR3454.txt fluxParams_HR9087.txt
do
echo "Running Cal on $FILE1 "
python testFluxCal.py $FILE1
done
