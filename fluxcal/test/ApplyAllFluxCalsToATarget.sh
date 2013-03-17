#! /bin/bash

#for FILE1 in fluxParams_G158-100.txt fluxParams_crab.txt fluxParams_crabNight2.txt fluxParams_nltt11748.txt fluxParams_nltt11748Night2.txt fluxParams_sdssj0926.txt
for FILE1 in fluxParams_sdssj0926.txt
do

#for FILE2 in fluxParams_hiltner600.txt fluxParams_HR3454.txt fluxParams_HR9087.txt fluxParams_feige66.txt fluxParams_bd25.txt fluxParams_hd849a.txt
for FILE2 in fluxParams_hd849a.txt
do
echo "Running Cal on $FILE1 with $FILE2"
python testApertureSpectrum.py $FILE1 $FILE2 0
done

#for FILE2 in fluxParams_HR3454_1.txt fluxParams_feige66_1.txt fluxParams_bd25_1.txt fluxParams_hd849a_1.txt
for FILE2 in fluxParams_hd849a_1.txt
do
echo "Running Cal on $FILE1 with $FILE2"
python testApertureSpectrum.py $FILE1 $FILE2 1
done

done

#for FILE1 in fluxParams_G158-100_1.txt fluxParams_crab_1.txt fluxParams_crabNight2_1.txt fluxParams_sdssj0926_1.txt
#do

#for FILE2 in fluxParams_hiltner600.txt fluxParams_HR3454.txt fluxParams_HR9087.txt fluxParams_feige66.txt fluxParams_bd25.txt fluxParams_hd849a.txt
#for FILE2 in fluxParams_hd849a.txt
#do
#echo "Running Cal on $FILE1 with $FILE2"
#python testApertureSpectrum.py $FILE1 $FILE2 2
#done

#for FILE2 in fluxParams_HR3454_1.txt fluxParams_feige66_1.txt fluxParams_bd25_1.txt fluxParams_hd849a_1.txt
#for FILE2 in fluxParams_hd849a_1.txt
#do
#echo "Running Cal on $FILE1 with $FILE2"
#python testApertureSpectrum.py $FILE1 $FILE2 3
#done

#done

