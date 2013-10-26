To use tsBinner, run this command in this directory:

python setup.py build_ext --inplace

and that will convert tsBinner.pyx to tsBinner.so after writing tsBinner.c


On a mac with enthought python, set this environment varible if the 
compilation has trouble finding numpy/arrayobject.h

export export CPATH=/Users/stoughto/Library/Enthought/Canopy_64bit/User/lib/python2.7/site-packages/numpy/core/include

where you need to replace the actual path with the actual path to the
directory that containd numpy/arrayobject.h
