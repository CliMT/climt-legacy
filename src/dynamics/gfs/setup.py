from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup (
        ext_modules = cythonize([
            Extension("_gfs_dynamics", ["_gfs_dynamics.pyx"],
               libraries=["gfsDycore","shtns_omp","lapack","fftw3_omp","fftw3","rt","m"] )
            ])

      )
