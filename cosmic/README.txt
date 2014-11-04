To use tsBinner, run this command in this directory:

python setup.py build_ext --inplace

and that will convert tsBinner.pyx to tsBinner.so after writing tsBinner.c


On a mac with enthought python, set this environment varible if the 
compilation has trouble finding numpy/arrayobject.h

export export CPATH=/Users/stoughto/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/numpy/core/include

where you need to replace the actual path with the actual path to the
directory that containd numpy/arrayobject.h


======

To run from a dict file, make a directory and go there.  
(For example, /Scratch/stoughto/geminga)

Run the cosmic finding like this, in the background, saving to a log file:

nohup  python -u $PYTHONPATH/cosmic/CosmicRunFromDictFile.py geminga-a.dict cosmic.dict > geminga-a.log 2>&1 &


The file geminga-a.dict defines the data set to use.  
The file cosmic.dict defins parameters to use in finding cosmic rays.

add a third parameter to run the "summary" and make plots and such.
