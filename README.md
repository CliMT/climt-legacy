# CliMT: Climate Modeling and Diagnostics Toolkit

CliMT is a Python-based software component toolkit which provides a
 flexible, multi-purpose problem-solving environment geared to climate
 science problems.

## Requirements:
* Python 2.4 or higher
* Numpy 1.0 or higher
* GNU gcc and a Fortran compiler

## Obtaining CliMT

Make sure you have a git client installed. Most operating systems
have a git client available, you can choose any.

To obtain CliMT, use the following command:

git clone --recursive https://github.com/CliMT/climt-python.git

Please note the **recursive** keyword. Without this, the dynamical
core will *not* be downloaded.

## Installation:
Make sure that numpy is installed, and that the version of f2py 
on your path is the one installed by numpy. Then:

**python setup.py build**

will build the package using the default fortran compiler. For
information on the compilers available on your platform:

**f2py -c --help-fcompiler**

To build with a different compiler:

**python setup.py build --fcompiler=<compiler>**

where <compiler> is one of the choices offered by f2py above.

When the build is complete:

**python setup.py install --prefix=<dir where you want climt to live>**

You can also install the "lite" version of CliMT, used in 
Ray Pierrehumbert's book "Principles of Planetary Climate",
by doing

**python setup.py install --lite**

Rodrigo Caballero (rodrigo.caballero@misu.su.se)
Feb 2012


