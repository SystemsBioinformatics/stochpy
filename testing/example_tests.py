#! /usr/bin/env python
"""
StochPy - Stochastic modeling in Python (http://stochpy.sourceforge.net)

Copyright (C) 2010-2025 T.R Maarlveld, B.G. Olivier, F.J. Bruggeman. All rights reserved.

This file:
Author: Brett G. Olivier (b.g.olivier@vu.nl)
Address: Vrije Universiteit Amsterdam, Amsterdam, Netherlands.


Permission to use, modify, and distribute this software is given under the
terms of the StochPy (BSD style) license. 

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
"""

print('# StochPy Basic tests')

print('## Import and load')

import stochpy
smod = stochpy.SSA()

print('\nDone.\n\n## Basic Simulation with the Direct method')

smod.DoStochSim(IsTrackPropensities=True)
smod.data_stochsim.simulation_endtime
smod.data_stochsim.simulation_timesteps
smod.GetWaitingtimes()
smod.PrintWaitingtimesMeans()

print('\nDone.\n\n## Do some Plotting')

smod.PlotSpeciesTimeSeries()
smod.PlotWaitingtimesDistributions()
smod.PlotPropensitiesTimeSeries()

print('\nDone.\n\n## Write data to a text file')

smod.Export2File()
smod.Export2File(analysis='distribution')
smod.Export2File(analysis='distribution',datatype='species')
smod.Export2File(analysis='mean',datatype='species')
smod.Export2File(analysis='std',datatype='species')
smod.Export2File(analysis='autocorrelation',datatype='species')

print('\nDone.\n\n## Show the means from the data of 3-th trajectory')

smod.DoStochSim(trajectories=3) # multiple trajectories
smod.data_stochsim.simulation_trajectory
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()

print('\nDone.\n\n## Switch to data from trajectory 1 and show the means of each species')

smod.GetTrajectoryData(1)
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()

print('\nDone.\n\n## Do one long simulation')

smod.DoStochSim(trajectories=1,end=1000000,mode='steps')
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()

print('\nDone.\n\n## Plot the PDF for different bin sizes')

smod.PlotSpeciesDistributions()
smod.PlotSpeciesDistributions(bin_size=5)  # larger bin size
smod.PlotSpeciesDistributions(bin_size=10) # again a larger bin size
smod.Export2File(analysis='distribution',datatype='species')

print('\nDone.\n\n## Usage of the Reload Function: `Ksyn = 20, kdeg = 0.2`')

smod.ChangeParameter('Ksyn',20.0)
smod.ChangeParameter('Kdeg',0.2)
smod.DoStochSim()
smod.PrintSpeciesMeans()   # should be ~Ksyn/Kdeg

print('\nDone.\n\n## Use another model to show the Interpolation features')

smod.Model('dsmts-001-01.xml.psc')
smod.DoStochSim(trajectories=1000,end=50,mode='time')
smod.GetRegularGrid(n_samples=51)
smod.PlotAverageSpeciesTimeSeries()
smod.PrintAverageSpeciesTimeSeries()
smod.Export2File(datatype='species',analysis='timeseries',IsAverage=True)

print('\nDone.\n\n## Test each method for different models:')

smod.Model('Autoreg.psc')
smod.DoStochSim(trajectories=1,end=1000,mode='steps')
smod.Method('NextReactionMethod')
smod.DoStochSim(trajectories=1,end=1000,mode='steps')
smod.data_stochsim.species
smod.PlotWaitingtimesDistributions()
smod.Method('FirstReactionMethod')
smod.DoStochSim(trajectories=1,end=1000,mode='steps')
smod.Method('TauLeaping')
smod.DoStochSim(trajectories=1,end=1000,mode='steps')

smod.Model('DecayingDimerizing.psc')
smod.DoStochSim(method = 'Direct',trajectories=1,end=50,mode='time')
smod.DoStochSim(method = 'NextReactionMethod',trajectories=1,end=50,mode='time')
smod.DoStochSim(method = 'FirstReactionMethod',trajectories=1,end=50,mode='time')
smod.PlotWaitingtimesDistributions()
smod.DoStochSim(method = 'TauLeaping',trajectories=1,end=50,mode='time',epsilon=0.03)  # Should outperform all other implementations
smod.PlotSpeciesTimeSeries()
#smod.PlotWaitingtimesDistributions()   # Should give an error

smod.Model('chain500.psc')
smod.DoStochSim(method = 'Direct',trajectories=1,end=10000,mode='steps')
smod.DoStochSim(method = 'NextReactionMethod',trajectories=1,end=10000,mode='steps') # should outperform the direct method and all other implementations

print('\n### Use the Next Reaction Method to test a model with a time event\n')

smod.Model('dsmts-003-03.xml.psc')
smod.DoStochSim(method = 'NextReactionMethod')
smod.DoTestsuite()

print('\n### Use the First Reaction method to test a model with a concentration event\n')

smod.Model('dsmts-003-04.xml.psc')
smod.DoStochSim(method = 'FirstReactionMethod')
smod.DoTestsuite()

print('\n### Volume Models\n')

smod.Model('dsmts-001-11.xml.psc')
smod.DoStochSim(method = 'Direct',trajectories=1000,end=50,mode ='time')
smod.PrintAverageSpeciesTimeSeries()
