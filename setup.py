#!/usr/bin/env python
from __future__ import division, print_function, absolute_import

__doc__ = "StochPy (Stochastic modeling in Python) provides several stochastic simulation algorithms to simulate (bio)chemical systems of reactions in a stochastic manner."
__version__ = "2.3"

import os

try:
    from numpy.distutils.core import setup
except Exception as ex:
    print(ex)
    print("StochPy requires NumPy\n")
    print("See http://numpy.scipy.org/ for more information about NumPy")
    os.sys.exit(-1)

local_path = os.path.dirname(os.path.abspath(os.sys.argv[0]))		# Get the dir of setup.py
os.chdir(local_path)

mydata_files = []
modfold = os.path.join(local_path, 'stochpy', 'pscmodels')
mods = os.listdir(modfold)

mypackages = ['stochpy','stochpy.lib','stochpy.modules','stochpy.pscmodels','stochpy.implementations','stochpy.core2','stochpy.tools']		# My subpackage list
mymodules = []

setup(name="StochPy",
    version = __version__,
    description = __doc__,
    long_description = """
    Welcome to the installation of StochPy {0:s}!

    StochPy (Stochastic modeling in Python) is a flexible software tool for stochastic simulation in cell biology. It provides various stochastic simulation algorithms, SBML support, analyses of the probability distributions of molecule copy numbers and event waiting times, analyses of stochastic time series, and a range of additional statistical functions and plotting facilities for stochastic simulations.
    """.format(__version__),
    author = "T.R. Maarleveld",
    author_email = "tmd200@users.sourceforge.net",
    maintainer = "T.R. Maarleveld",
    maintainer_email = "tmd200@users.sourceforge.net",
    url = "http://stochpy.sourceforge.net",
    download_url = "http://stochpy.sourceforge.net/download.html",
    license = " BSD License ",
    keywords = " Bioinformatics, Computational Systems Biology, Bioinformatics, Modeling, Simulation, Stochastic Simulation Algorithms, Stochastic", 
    zip_safe = False,
    requires = ['NumPy'],
    platforms = ["Windows", "Linux","Mac OS-X"],#, "Solaris", "", "Unix"],
    classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: End Users/Desktop',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: BSD License',
    'Natural Language :: English',
    'Operating System :: OS Independent',    
    'Programming Language :: Python :: 2.6',
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3.4',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Topic :: Scientific/Engineering :: Simulations'],
    packages = mypackages,
    data_files = mydata_files
    )
