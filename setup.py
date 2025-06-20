#!/usr/bin/env python
from __future__ import division, print_function, absolute_import
__doc__ = "StochPy (Stochastic modeling in Python) provides several stochastic simulation algorithms to simulate (bio)chemical systems of reactions in a stochastic manner."
__version__ = "2.5"
import os
import sys
import time

local_path = os.path.dirname(os.path.abspath(os.sys.argv[0]))


from setuptools import setup, find_packages, find_namespace_packages

install_requires_src = ['numpy', 'packaging', 'python_libsbml', 'lxml', 'scipy', 'matplotlib', 'ipython']
extras_require_src = {'all': ['numpy', 'packaging', 'python_libsbml', 'lxml', 'scipy', 'matplotlib', 'ipython']}
tests_require_src = ['numpy', 'packaging', 'python_libsbml', 'lxml', 'scipy', 'matplotlib', 'ipython']
mydata_files = []

#mypackages = ['stochpy','stochpy.lib','stochpy.modules','stochpy.pscmodels','stochpy.implementations','stochpy.core2','stochpy.tools']  # My subpackage list
mypackages = find_packages(where='stochpy',
                           include=['lib', 'modules', 'pscmodels', 'implementations', 'core2', 'tools'])

setup(
    package_dir={"stochpy": "stochpy"},
    packages=mypackages,
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
    install_requires=install_requires_src,
    extras_require=extras_require_src,
    python_requires='>=3.9',
    platforms=["Windows", "Linux", "Mac"],  # "Solaris", "", "Unix"],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Development Status :: 6 - Mature',
        'Environment :: Console',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    data_files=mydata_files
)
