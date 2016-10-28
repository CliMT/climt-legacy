from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran convect43c.f90 -c -o emanuel_fortran.o -O3 -fPIC'
print fortran_mod_comp
system(fortran_mod_comp)

ext_modules = [Extension(# module name:
                         '_emanuel_convection_new',
                         # source file:
                         ['_emanuel_convection.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3', '-lgfortran'],
                         # other files to link to
                         extra_link_args=['emanuel_fortran.o', '-lgfortran'])]

setup(name = '_emanuel_convection_new',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)

