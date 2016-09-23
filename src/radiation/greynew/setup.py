from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking

ext_modules = [Extension(# module name:
                         '_new_grey',
                         # source file:
                         ['_grey_gas.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3', '-lgfortran'],
                         # other files to link to
                         )]

setup(name = '_new_grey',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)
