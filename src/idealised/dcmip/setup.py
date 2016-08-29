from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
# This line only needed if building with NumPy in Cython file.
from numpy import get_include
from os import system

# compile the fortran modules without linking
fortran_mod_comp = 'gfortran dcmip_initial_conditions_test_4_v3.f90 -c -o dcmip_test4.o -O3 -fPIC'
print fortran_mod_comp
system(fortran_mod_comp)
fortran_mod_comp = 'gfortran dcmip_initial_conditions_test_5_v1.f90 -c -o dcmip_test5.o -O3 -fPIC'
print fortran_mod_comp
system(fortran_mod_comp)

ext_modules = [Extension(# module name:
                         'dcmip',
                         # source file:
                         ['dcmip_initial_conditions.pyx'],
                         # other compile args for gcc
                         extra_compile_args=['-fPIC', '-O3', '-lgfortran'],
                         # other files to link to
                         extra_link_args=['dcmip_test4.o','dcmip_test5.o', '-lgfortran'])]

setup(name = 'dcmip',
      cmdclass = {'build_ext': build_ext},
      # Needed if building with NumPy.
      # This includes the NumPy headers when compiling.
      include_dirs = [get_include()],
      ext_modules = ext_modules)
