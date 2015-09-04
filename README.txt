StochPy Stochastic modeling in Python
=====================================

Copyright (c) 2011-2015, Timo R. Maarleveld, Brett G. Olivier, and Frank J. Bruggeman
All rights reserved.

StochPy is distributed under a BSD style licence.

Author information
------------------

Timo R. Maarleveld, Brett G. Olivier, and Frank J. Bruggeman
Centrum Wiskunde en Informatica, Amsterdam, Netherlands
VU University, Amsterdam, Netherlands

e-mail: tmd200@users.sourceforge.net
web: http://sourceforge.net/projects/stochpy/

Documentation can be found in the user guide (see Documentation directory or http://stochpy.sourceforge.net/html/userguide.html) 

Installation
------------
The following software is required before installling StochPy (see user guide for more details):

- Python 2.6+ or Python 3.4+
- NumPy 1.x+
- Matplotlib (optional)
- libsbml (optional)
- libxml2 (optional)
- mpmath (optional)

Linux/MAC OS/Cygwin
~~~~~~~~~~~~~~~~~~~

1) cd to directory StochPy-2.1.0
2) sudo python setup.py install

Windows
~~~~~~~
Use the available windows installer or the setup file

Usage
-----

import stochpy
smod = stochpy.SSA()

# 1: Basic Simulation with the Direct method
smod.DoStochSim(IsTrackPropensities=True)
smod.data_stochsim.simulation_endtime
smod.data_stochsim.simulation_timesteps
smod.GetWaitingtimes()
smod.PrintWaitingtimesMeans()
# 2: Do some Plotting
smod.PlotSpeciesTimeSeries()
smod.PlotWaitingtimesDistributions()
smod.PlotPropensitiesTimeSeries()
# 3: Write data to a text file
smod.Export2File()
smod.Export2File(analysis='distribution')
smod.Export2File(analysis='distribution',datatype='species')
smod.Export2File(analysis='mean',datatype='species')
smod.Export2File(analysis='std',datatype='species')
smod.Export2File(analysis='autocorrelation',datatype='species')
# 4: Show the means from the data of 3-th trajectory
smod.DoStochSim(trajectories=3) # multiple trajectories
smod.data_stochsim.simulation_trajectory
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()
# 5: Switch to data from trajectory 1 and show the means of each species
smod.GetTrajectoryData(1)
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()
# 6: Do one long simulation
smod.DoStochSim(trajectories=1,end=1000000,mode='steps')
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()
# 7: Plot the PDF for different bin sizes
smod.PlotSpeciesDistributions()
smod.PlotSpeciesDistributions(bin_size=5)  # larger bin size
smod.PlotSpeciesDistributions(bin_size=10) # again a larger bin size
smod.Export2File(analysis='distribution',datatype='species')

# 8: Usage of the Reload Function: Ksyn = 20, kdeg = 0.2
smod.ChangeParameter('Ksyn',20.0)
smod.ChangeParameter('Kdeg',0.2)
smod.DoStochSim()
smod.PrintSpeciesMeans()   # should be ~Ksyn/Kdeg

# 9: Use another model to show the Interpolation features
smod.Model('dsmts-001-01.xml.psc')
smod.DoStochSim(trajectories=1000,end=50,mode='time') 
smod.GetRegularGrid(npoints=51)
smod.PlotAverageSpeciesTimeSeries()
smod.PrintAverageSpeciesTimeSeries()
smod.Export2File(datatype='species',analysis='timeseries',IsAverage=True)

# 9: Test each method for different models: 
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

# 10: Use the Next Reaction Method to test a model with a time event
smod.Model('dsmts-003-03.xml.psc') 
smod.DoStochSim(method = 'NextReactionMethod')
smod.DoTestsuite()

# 11: Use the First Reaction method to test a model with a concentration event 
smod.Model('dsmts-003-04.xml.psc')
smod.DoStochSim(method = 'FirstReactionMethod')
smod.DoTestsuite()

# 12: Volume Models
smod.Model('dsmts-001-11.xml.psc') 
smod.DoStochSim(method = 'Direct',trajectories=1000,end=50,mode ='time')
smod.PrintAverageSpeciesTimeSeries()
