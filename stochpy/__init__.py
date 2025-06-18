#! /usr/bin/env python
"""
StochPy - Stochastic modeling in Python (http://stochpy.sourceforge.net)

Copyright (C) 2010-2015 T.R Maarlveld, B.G. Olivier, F.J. Bruggeman all rights reserved.

Timo R. Maarleveld (tmd200@users.sourceforge.net)
Centrum Wiskunde & Informatica, Amsterdam, Netherlands
VU University, Amsterdam, Netherlands

Permission to use, modify, and distribute this software is given under the
terms of the StochPy (BSD style) license.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
"""
from __future__ import division, print_function, absolute_import

__doc__ =   """
            StochPy: Stochastic Modeling in Python
            =====================================

            StochPy (Stochastic modeling in Python) is a flexible software tool for stochastic simulation in cell biology. It provides various stochastic
            simulation algorithms, SBML support, analyses of the probability distributions of molecule copy numbers and event waiting times, analyses of stochastic time
            series, and a range of additional statistical functions and plotting facilities for stochastic simulations.

            Options:
            --------
            - Stochastic Simulations
            - Variety of stochastic simulation output analysis:
              --> Time Simulation
              --> Distribution
              --> Waiting times
              --> Propensities
            - Cell Division simulations
            - SBML and PSC MDL input format.

            StochPy can be used in an interactive Python shell:

            Usage
            -----
            >>> import stochpy
            >>> utils = stochpy.Utils()
            >>> utils.doExample1()
            >>> utils.doExample2()
            >>> smod = stochpy.SSA()   # stochastic simulation algorithm module
            >>> help(smod)
            >>> help(stochpy.SSA)      # (some windows versions)
            >>> stochpy?
            >>> smod.DoStochSim()
            >>> smod.PlotSpeciesTimeSeries()
            >>> converter = stochpy.SBML2PSC()
            >>> converter??
            >>> help(stochpy.SBML2PSC)
            """

from .core2.version import __version__

import os, shutil, sys

try:
    import readline
    _IsReadline = True
except:
    _IsReadline = False

try:
    # from numpy.distutils.core import setup, Extension
    import numpy
    _IsNumPy = True
except Exception as ex:
    _IsNumPy = False
    print(ex)
    print("StochPy requires NumPy")
    print("See http://numpy.scipy.org/ for more information about NumPy")
    os.sys.exit(-1)

try:
    import matplotlib
    _IsMPL = True
except:
    _IsMPL = False
    print("Warning: The Matplotlib module is not available, so plotting is not possible")
    print("Info: See http://matplotlib.sourceforge.net/ for more information about Matplotlib.")

_IsPlotting = False
try:
    import matplotlib.pyplot as plt
    _IsPlotting = True
except Exception as er:
    print(er)

def InitiateModels(directory):
    """
    Build several models written in PSC MDL and SBML

    Input:
     - *directory* (string)
    """
    from .pscmodels import Burstmodel
    from .pscmodels import BirthDeath
    from .pscmodels import ImmigrationDeath
    from .pscmodels import DecayingDimerizing
    from .pscmodels import Autoreg
    from .pscmodels import CellDivision as celldivision
    from .pscmodels import GeneDuplication
    from .pscmodels import dsmts_001_01
    from .pscmodels import dsmts_001_11
    from .pscmodels import dsmts_001_19
    from .pscmodels import dsmts_002_10
    from .pscmodels import dsmts_003_03
    from .pscmodels import dsmts_003_04
    from .pscmodels import chain5
    from .pscmodels import chain50
    from .pscmodels import chain500
    from .pscmodels import chain1500
    from .pscmodels import Isomerization
    from .pscmodels import Polymerase
    from .pscmodels import TranscriptionIntermediate
    from .pscmodels import Schlogl
    from .pscmodels import SignalingTimeVaryingL
    from .pscmodels import Signaling3cCD

    models = {}
    models['Signaling3cCD.psc'] = Signaling3cCD.model
    models['SignalingTimeVaryingL.psc'] = SignalingTimeVaryingL.model
    models['Schlogl.psc'] = Schlogl.model
    models['Burstmodel.psc'] = Burstmodel.model
    models['ImmigrationDeath.psc'] = ImmigrationDeath.model
    models['BirthDeath.psc'] = BirthDeath.model
    models['DecayingDimerizing.psc'] = DecayingDimerizing.model
    models['Autoreg.psc'] = Autoreg.model
    models['Autoreg.xml'] = Autoreg.xml_model
    models['CellDivision.psc'] = celldivision.model
    models['GeneDuplication.psc'] = GeneDuplication.model
    models['Isomerization.psc'] = Isomerization.model
    models['Polymerase.psc'] = Polymerase.model
    models['TranscriptionIntermediate.psc'] = TranscriptionIntermediate.model
    models['dsmts-001-01.xml.psc'] = dsmts_001_01.model
    models['dsmts-001-01.xml'] = dsmts_001_01.xml_model
    models['dsmts-001-11.xml.psc'] = dsmts_001_11.model
    models['dsmts-001-11.xml'] = dsmts_001_11.xml_model
    models['dsmts-001-19.xml.psc'] = dsmts_001_19.model
    models['dsmts-001-19.xml'] = dsmts_001_19.xml_model
    models['dsmts-002-10.xml.psc'] = dsmts_002_10.model
    models['dsmts-002-10.xml'] = dsmts_002_10.xml_model
    models['dsmts-003-03.xml.psc'] = dsmts_003_03.model
    models['dsmts-003-03.xml'] = dsmts_003_03.xml_model
    models['dsmts-003-04.xml.psc'] = dsmts_003_04.model
    models['dsmts-003-04.xml'] = dsmts_003_04.xml_model
    models['chain5.psc'] = chain5.model
    models['chain50.psc'] = chain50.model
    models['chain500.psc'] = chain500.model
    models['chain1500.psc'] = chain1500.model

    model_names = list(models)
    dir_models = os.listdir(directory)
    for mod_name in model_names:
        if mod_name not in dir_models:
            print("Info: Model {0:s} copied to {1:s}".format(mod_name ,directory) )
            file_out = open(os.path.join(directory,mod_name),'w')
            file_out.write(models[mod_name])
            file_out.close()


output_dir = None
model_dir = None
if os.sys.platform != 'win32':
    if not os.path.exists(os.path.join(os.path.expanduser('~'),'Stochpy')):
        os.makedirs(os.path.join(os.path.expanduser('~'),'Stochpy'))
    if not os.path.exists(os.path.join(os.path.expanduser('~'),'Stochpy', 'pscmodels')):
        os.makedirs(os.path.join(os.path.expanduser('~'),'Stochpy','pscmodels'))
    if not os.path.exists(os.path.join(os.path.expanduser('~'),'Stochpy', 'temp')):
        os.makedirs(os.path.join(os.path.expanduser('~'),'Stochpy','temp'))

    output_dir = os.path.join(os.path.expanduser('~'),'Stochpy')
    model_dir = os.path.join(os.path.expanduser('~'),'Stochpy','pscmodels')
    temp_dir = os.path.join(os.path.expanduser('~'),'Stochpy','temp')
    InitiateModels(model_dir)
else:
    if not os.path.exists(os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy')):
        os.makedirs(os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy'))
    if not os.path.exists(os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy','pscmodels')):
        os.makedirs(os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy','pscmodels'))
    if not os.path.exists(os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy','temp')):
        os.makedirs(os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy','temp'))

    output_dir = os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy',)
    model_dir = os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy','pscmodels')
    temp_dir = os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy','temp')
    InitiateModels(model_dir)

from .modules.SBML2PSC import SBML2PSC
from .modules.StochSim import SSA
from .modules.StochPyUtils import Utils
from .modules.StochPyCellDivision import CellDivision
from .modules.StochPyDemo import Demo
from .modules import Analysis as Analysis

try:
    from .modules.NucleosomeTool import NucleosomeModelBuilder
    from .modules.NucleosomeTool import NucleosomeSimulator
except Exception as er:
    pass # ignore

def DeletePreviousOutput(path,type):
    """
    Delete output of earlier simulations

    Input:
     - *path* (string)
     - *type* (string)
    """
    for filename in os.listdir(path):
        if filename.endswith(type):
            filename_path = os.path.join(path,filename)
            os.remove(filename_path)

def DeleteExistingData(path):
    """
    Delete all existing StochKit simulation data

    Input:
     - *path* (string)
    """
    if os.path.exists(path):
        for maps in os.listdir(path):
            dir2delete = os.path.join(path,maps)
            shutil.rmtree(dir2delete, ignore_errors=True)

def SaveInteractiveSession(filename='interactiveSession.py',path=output_dir):
    """
    Save the interactive session

    Input:
     - *filename*: [default = interactiveSession.py'] (string)
     - *path*: (string)
    """
    if not _IsReadline:
        print("Error: install 'readline' first")
    elif _IsReadline:
        historyPath = os.path.join(path,filename)
        if not os.path.exists(path):
            os.makedirs(path)

        readline.write_history_file(historyPath)
        file_in = open(historyPath,'r')
        history_list = file_in.readlines()
        n_import_statement = 0
        for command in history_list:
            if 'import' in command and 'stochpy' in command:
                n_import_statement += 1

        n=0
        file_out = open(historyPath,'w')
        for command in history_list:
            if 'import' in command and 'stochpy' in command:
                n+=1
            if n==n_import_statement:
                file_out.write(command)
        file_out.close()
        print("Info: Interactive session successfully saved at {0:s}".format(historyPath) )
        print("Info: use 'ipython {0:s} to restart modeling with this interactive session".format(filename) )


DeletePreviousOutput(temp_dir,'.dat')
DeletePreviousOutput(temp_dir,'.xml')
DeletePreviousOutput(temp_dir,'.txt')
DeletePreviousOutput(temp_dir,'temp_parse_module')
DeleteExistingData(temp_dir)
#readline.clear_history()

print("""
#######################################################################
#                                                                     #
#            Welcome to the interactive StochPy environment           #
#                                                                     #
#######################################################################
#  StochPy: Stochastic modeling in Python                             #
#  http://stochpy.sourceforge.net                                     #
#  Copyright(C) T.R Maarleveld, B.G. Olivier, F.J Bruggeman 2010-2015 #
#  DOI: 10.1371/journal.pone.0079345                                  #
#  Email: tmd200@users.sourceforge.net                                #
#  VU University, Amsterdam, Netherlands                              #
#  Centrum Wiskunde Informatica, Amsterdam, Netherlands               #
#  StochPy is distributed under the BSD licence.                      #
#######################################################################
""")
print("Version {0:s}".format(__version__) )
print("Output Directory: {0:s}".format(output_dir) )
print("Model Directory: {0:s}".format(model_dir) )
#print("Warning: Figure freezing? Try a different matplotlib backend (stochpy.plt.switch_backend) and/or set IsInteractive to False (see user guide)")
