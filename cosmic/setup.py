from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("tsBinner",["tsBinner.pyx"])]

setup(
  name ='time stamp binner',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
