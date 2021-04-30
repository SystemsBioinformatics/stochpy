# StochPy Stochastic modeling in Python

StochPy is a versatile stochastic modeling package which is designed for stochastic simulation of molecular control networks

* File releases: http://sourceforge.net/projects/stochpy
* Source code: https://github.com/SystemsBioinformatics/stochpy

StochPy is open source software distributed under the BSD 3-Clause License, see LICENSE file for more details.

## Documentation

Documentation can be found in the user guide (see Documentation directory or in [sourceforge](http://stochpy.sourceforge.net/html/userguide.html))

## Installation

The following software is required before installing StochPy (see user guide for more details):

- Python 2.6+ or Python 3.4+
- [NumPy 1.x+](http://www.numpy.org)
- [Matplotlib](https://matplotlib.org) (optional)
- [libsbml](http://sbml.org/Software/libSBML) (optional)
- [libxml2](http://xmlsoft.org) (optional)
- [mpmath](http://mpmath.org) (optional)

Install the StochPy dependencies with PIP using the following command (in your StochPy Python environment):
```bash
pip install numpy matplotlib python-libsbml mpmath lxml
```

#### Linux/MAC OS/Cygwin


In the directory where you downloaded StochPy, go to the directory _StochPy-2.1.0_ and exec the next command:

```bash
sudo python setup.py install
```

### Windows

Use the available windows installer or the setup file

## Getting Started

You can run `ipython` and import `stochpy`:

```py
import stochpy
smod = stochpy.SSA()
```

### Basic Simulation with the Direct method
```py
smod.DoStochSim(IsTrackPropensities=True)
smod.data_stochsim.simulation_endtime
smod.data_stochsim.simulation_timesteps
smod.GetWaitingtimes()
smod.PrintWaitingtimesMeans()
```

### Do some Plotting
```py
smod.PlotSpeciesTimeSeries()
smod.PlotWaitingtimesDistributions()
smod.PlotPropensitiesTimeSeries()
```

### Write data to a text file
```py
smod.Export2File()
smod.Export2File(analysis='distribution')
smod.Export2File(analysis='distribution',datatype='species')
smod.Export2File(analysis='mean',datatype='species')
smod.Export2File(analysis='std',datatype='species')
smod.Export2File(analysis='autocorrelation',datatype='species')
```

### Show the means from the data of 3-th trajectory
```py
smod.DoStochSim(trajectories=3) # multiple trajectories
smod.data_stochsim.simulation_trajectory
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()
```

### Switch to data from trajectory 1 and show the means of each species
```py
smod.GetTrajectoryData(1)
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()
```

### Do one long simulation
```py
smod.DoStochSim(trajectories=1,end=1000000,mode='steps')
smod.PrintSpeciesMeans()
smod.PrintSpeciesStandardDeviations()
```

### Plot the PDF for different bin sizes
```py
smod.PlotSpeciesDistributions()
smod.PlotSpeciesDistributions(bin_size=5)  # larger bin size
smod.PlotSpeciesDistributions(bin_size=10) # again a larger bin size
smod.Export2File(analysis='distribution',datatype='species')
```

### Usage of the Reload Function: `Ksyn = 20, kdeg = 0.2`
```py
smod.ChangeParameter('Ksyn',20.0)
smod.ChangeParameter('Kdeg',0.2)
smod.DoStochSim()
smod.PrintSpeciesMeans()   # should be ~Ksyn/Kdeg
```

### Use another model to show the Interpolation features
```py
smod.Model('dsmts-001-01.xml.psc')
smod.DoStochSim(trajectories=1000,end=50,mode='time')
smod.GetRegularGrid(npoints=51)
smod.PlotAverageSpeciesTimeSeries()
smod.PrintAverageSpeciesTimeSeries()
smod.Export2File(datatype='species',analysis='timeseries',IsAverage=True)
```

### Test each method for different models:
```py
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
```
```py
smod.Model('DecayingDimerizing.psc')
smod.DoStochSim(method = 'Direct',trajectories=1,end=50,mode='time')
smod.DoStochSim(method = 'NextReactionMethod',trajectories=1,end=50,mode='time')
smod.DoStochSim(method = 'FirstReactionMethod',trajectories=1,end=50,mode='time')
smod.PlotWaitingtimesDistributions()
smod.DoStochSim(method = 'TauLeaping',trajectories=1,end=50,mode='time',epsilon=0.03)  # Should outperform all other implementations
smod.PlotSpeciesTimeSeries()
#smod.PlotWaitingtimesDistributions()   # Should give an error
```
```py
smod.Model('chain500.psc')
smod.DoStochSim(method = 'Direct',trajectories=1,end=10000,mode='steps')
smod.DoStochSim(method = 'NextReactionMethod',trajectories=1,end=10000,mode='steps') # should outperform the direct method and all other implementations
```

### Use the Next Reaction Method to test a model with a time event
```py
smod.Model('dsmts-003-03.xml.psc')
smod.DoStochSim(method = 'NextReactionMethod')
smod.DoTestsuite()
```

### Use the First Reaction method to test a model with a concentration event
```py
smod.Model('dsmts-003-04.xml.psc')
smod.DoStochSim(method = 'FirstReactionMethod')
smod.DoTestsuite()
```

### Volume Models
```py
smod.Model('dsmts-001-11.xml.psc')
smod.DoStochSim(method = 'Direct',trajectories=1000,end=50,mode ='time')
smod.PrintAverageSpeciesTimeSeries()
```

## Author information


Timo R. Maarleveld, Brett G. Olivier, and Frank J. Bruggeman
Centrum Wiskunde en Informatica, Amsterdam, Netherlands
VU University, Amsterdam, Netherlands

> e-mail: tmd200@users.sourceforge.net

## Publication

StochPy: A Comprehensive, User-Friendly Tool for Simulating Stochastic Biological Processes
http://dx.doi.org/10.1371/journal.pone.0079345

## Licence
Copyright (c) 2011-2021, Timo R. Maarleveld, Brett G. Olivier, and Frank J. Bruggeman
Vrije Universiteit Amsterdam. All rights reserved.

StochPy is open source software distributed under the BSD 3-Clause License see LICENSE file for more details.