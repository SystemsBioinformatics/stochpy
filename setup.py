#!/usr/bin/env python
from __future__ import division, print_function, absolute_import
__doc__ = "StochPy (Stochastic modeling in Python) provides several stochastic simulation algorithms to simulate (bio)chemical systems of reactions in a stochastic manner."
__version__ = "2.4"
import os
import sys

try:
    import setuptools
    print('Using setuptools.')
except:
    print('Using distutils')

# not needed anympre
"""
try:
    import numpy
except Exception as ex:
    print(ex)
    print("StochPy requires NumPy\n")
    print("See http://numpy.scipy.org/ for more information about NumPy")
    sys.exit(-1)
"""

# from numpy.distutils.core import setup
from setuptools import setup

local_path = os.path.dirname(os.path.abspath(sys.argv[0]))  # Get the dir of setup.py
os.chdir(local_path)

# modfold = os.path.join(local_path, 'stochpy', 'pscmodels')
# mods = os.listdir(modfold)

# we now leave this to the pyproject definition
# mypackages = ['stochpy','stochpy.lib','stochpy.modules','stochpy.pscmodels','stochpy.implementations','stochpy.core2','stochpy.tools']  # My subpackage list
mydata_files = []
mymodules = []
mypackages = []

setup(
    name="StochPy",
    version=__version__,
    description=__doc__,
    long_description="""
    Welcome to the installation of StochPy {0:s}!
    StochPy (Stochastic modeling in Python) is a flexible software tool for stochastic simulation in cell biology. It provides various stochastic simulation algorithms, SBML support, analyses of the probability distributions of molecule copy numbers and event waiting times, analyses of stochastic time series, and a range of additional statistical functions and plotting facilities for stochastic simulations.
    """.format(__version__),
    author="T.R. Maarleveld",
    author_email="tmd200@users.sourceforge.net",
    maintainer="Brett Olivier",
    maintainer_email="tmd200@users.sourceforge.net",
    url="http://stochpy.sourceforge.net",
    download_url="https://github.com/SystemsBioinformatics/stochpy",
    license="BSD License",
    keywords="Bioinformatics, Computational Systems Biology, Bioinformatics, Modeling, Simulation, Stochastic Simulation Algorithms, Stochastic",
    zip_safe=False,
    python_requires='>=3.7',
    install_requires=[
        'numpy>=1.21.6',
        'scipy',
        'matplotlib',
        'lxml',
        'mpmath',
        'ipython',
    ],
    platforms=["Windows", "Linux", "Mac OS-X"],  # "Solaris", "", "Unix"],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Development Status :: 6 - Mature',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    packages=mypackages,
    data_files=mydata_files
)
