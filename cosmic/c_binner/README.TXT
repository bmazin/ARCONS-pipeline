For this statement to work from python:

from cosmic import c_binner

you need to do the following in this directory, after checking out from git

cd build
cmake ..
make

This will make the files c_binner.py and _c_binner.so in the directory
cosmic/c_binner/build

The __init__.py files in cosmic and cosmic/c_binner make this happen.

