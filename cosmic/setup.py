from distutils.core import setup
#from numpy.distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_modules = [
    Extension("tsBinner",
              ["tsBinner.pyx"],
              include_dirs=[np.get_include()]),
    ]

setup(
  name ='time stamp binner',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
