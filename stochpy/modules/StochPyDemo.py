 #! /usr/bin/env python
"""
Demo of StochPy Functionalities
===============================

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: August 06, 2015
"""
from __future__ import division, print_function, absolute_import

import stochpy
smod = stochpy.SSA()

try: input = raw_input # raw_input is renamed to input in python 3.x
except NameError: pass

class Demo():
    def __init__(self,DoSimulations = True):  
        """
        Demo class, by default all demo's are simulated.
        
        Input:
         - *DoSimulations* (Boolean)
        """
        input("press any button to continue\n")   
        if DoSimulations:
            self.DoDemoSimulations()
       
    def Demo1(self):
        """ Use the Immigration-Death model for doing basic simulations with the direct method, some plotting and exportation of results """    
        print("\n### (1) Basic simulations with the Immigration-Death model ###")
        input("press any button to continue\n")
        print(">>> smod = stochpy.SSA() # start SSA module")
        print(">>> smod.DoStochSim(IsTrackPropensities=True)")
        smod.DoStochSim(IsTrackPropensities=True)
        print(">>> smod.PrintWaitingtimesMeans()")
        smod.PrintWaitingtimesMeans()          
        print(">>> smod.PlotSpeciesTimeSeries() # plot time series of species")
        smod.PlotSpeciesTimeSeries()
        print(">>> smod.PlotPropensitiesTimeSeries() # plot time series of propensities")
        smod.PlotPropensitiesTimeSeries()
        print(">>> smod.PlotWaitingtimesDistributions()")
        smod.PlotWaitingtimesDistributions()
        print(">>> smod.Export2File()")
        smod.Export2File()
        print(">>> smod.Export2File(analysis='timeseries',datatype='propensities')")
        smod.Export2File(analysis='timeseries',datatype='propensities')
        print(">>> smod.Export2File(analysis='distribution',datatype='waitingtimes')")
        smod.Export2File(analysis='distribution',datatype='waitingtimes')
        
    def Demo2(self):    
        """ Use the Immigration-Death model for a stochastic simulation with multiple trajectories """
        print("\n### (2) Multiple trajectories ###")
        input("press any button to continue\n") 
        print(">>> smod.DoStochSim(trajectories=3,end=3000,mode='steps') # multiple trajectories")
        smod.DoStochSim(trajectories=3,end=3000) 
        print(">>> smod.PlotSpeciesTimeSeries()")
        smod.PlotSpeciesTimeSeries()
        print(">>> smod.data_stochsim.simulation_trajectory # trajectory number")
        smod.data_stochsim.simulation_trajectory
        print(">>> smod.PrintSpeciesMeans()")
        smod.PrintSpeciesMeans()
        print(">>> smod.PrintSpeciesStandardDeviations")
        smod.PrintSpeciesStandardDeviations()  
        print(">>> smod.GetTrajectoryData(1) # switch to data from trajectory 1")
        smod.GetTrajectoryData(1)      
        print(">>> smod.data_stochsim.simulation_trajectory # trajectory number")
        smod.data_stochsim.simulation_trajectory
        print(">>> smod.PrintSpeciesMeans()")
        smod.PrintSpeciesMeans()
        print(">>> smod.PrintSpeciesStandardDeviations")
        smod.PrintSpeciesStandardDeviations()  
    
    def Demo3(self):
        """ Use the Immigration-Death model to demonstrate probability density functions """
        print("\n### (3) Probability density functions (100000 time points) ###")
        input("press any button to continue\n")     
        print(">>> smod.DoStochSim(trajectories=1,end=1000000,mode='steps',IsTrackPropensities=1)")
        smod.DoStochSim(trajectories=1,end=100000,mode='steps',IsTrackPropensities=1)       
        print(">>> smod.PlotSpeciesDistributions() # plot species distributions")
        smod.PlotSpeciesDistributions()
        print(">>> smod.PlotPropensitiesDistributions('R2')")
        smod.PlotPropensitiesDistributions('R2')
        
    def Demo4(self):
        """ Use Birth-Death model to illustrate averaging multiple simulations """
        print("\n### (4) Averaging multiple simulations with the Birth-Death model ###")
        input("press any button to continue\n")     
        print(">>> smod.Model('dsmts-001-01.xml.psc') # parse a different model")
        smod.Model('dsmts-001-01.xml.psc')
        print(">>> smod.DoStochSim(trajectories=1000,end=50,mode='time')")
        smod.DoStochSim(trajectories=1000,end=50,mode='time') 
        print(">>> smod.GetRegularGrid(n_samples=51)")
        smod.GetRegularGrid(n_samples=51)
        print(">>> smod.PrintAverageSpeciesTimeSeries()")
        smod.PrintAverageSpeciesTimeSeries()  
        print(">>> smod.PlotAverageSpeciesTimeSeries()")
        smod.PlotAverageSpeciesTimeSeries()
        print(">>> smod.Export2File(analysis='timeseries',datatype='species', IsAverage = True)")
        smod.Export2File(analysis='timeseries',datatype='species', IsAverage = True)
        
    def Demo5(self):
        """ Use Decaying-Dimerizing model to illustrate performance differences between different stochastic simulation algorithms """
        print("\n### (5) Use Decaying-Dimerizing model to illustrate performance differences between different stochastic simulation algorithms ###")
        input("press any button to continue\n")
        print(">>> smod.Model('DecayingDimerizing.psc')")
        smod.Model('DecayingDimerizing.psc')
        print(">>> smod.DoStochSim(method = 'direct',trajectories=1,end=50,mode='time')")
        smod.DoStochSim(method = 'direct',trajectories=1,end=50,mode='time')
        print(">>> smod.PlotWaitingtimesDistributions()")
        smod.PlotWaitingtimesDistributions()     
        print(">>> smod.DoStochSim(method = 'Tauleap',trajectories=1,end=50,mode='time',epsilon=0.03) # should outperform all other implementations")
        smod.DoStochSim(method = 'Tauleap',trajectories=1,end=50,mode='time',epsilon=0.03)
        print(">>> smod.PlotSpeciesTimeSeries() # plot time series of species")
        smod.PlotSpeciesTimeSeries()     
        print(">>> stochpy.plt.xscale('log')")
        stochpy.plt.xscale('log')          
         
    def Demo6(self):
        """ Demo of StochPy's next reaction method handling events"""
        print("\n### (6) Demo of StochPy's next reaction method handling events ###")
        input("press any button to continue\n")
        print(">>> smod.Model('dsmts-003-03.xml.psc')")
        smod.Model('dsmts-003-03.xml.psc') 
        print(">>> smod.DoStochSim(method = 'NextReactionMethod',trajectories=1000,end=50,mode='time')")
        smod.DoStochSim(method = 'NextReactionMethod',trajectories=1000,end=50,mode='time')    
        print(">>> smod.GetRegularGrid()")
        smod.GetRegularGrid()
        print(">>> smod.PlotAverageSpeciesTimeSeries()")        
        smod.PlotAverageSpeciesTimeSeries()

    def Demo7(self):
        """ Demo of StochPy's first reaction method handling events"""
        print("\n### (7) Demo of StochPy's first reaction method handling events ###")
        input("press any button to continue\n")
        print(">>> smod.Model('dsmts-003-04.xml.psc')")
        smod.Model('dsmts-003-04.xml.psc') 
        print(">>> smod.DoStochSim(method = 'FirstReactionMethod',trajectories=1000,end=50,mode='time')")
        smod.DoStochSim(method = 'FirstReactionMethod',trajectories=1000,end=50,mode='time')    
        print(">>> smod.GetRegularGrid()")
        smod.GetRegularGrid()
        print(">>> smod.PlotAverageSpeciesTimeSeries()")        
        smod.PlotAverageSpeciesTimeSeries()

    def Demo8(self):
        """ Demo of StochPy's direct method for supporting events """
        print("\n### (8) Demo of StochPy's direct method handling events ###")
        input("press any button to continue\n")
        print(">>> smod.Model('dsmts-02-10.xml.psc')")
        smod.Model('dsmts-002-10.xml.psc') 
        print(">>> smod.DoStochSim(method = 'direct',trajectories=1000,end=50,mode='time')")
        smod.DoStochSim(method = 'direct',trajectories=1000,end=50,mode='time')    
        print(">>> smod.GetRegularGrid()")
        smod.GetRegularGrid()        
        print(">>> smod.PlotAverageSpeciesTimeSeries()")
        smod.PlotAverageSpeciesTimeSeries()   
    
    def Demo9(self): 
        """  Demo of StochPy's direct method for handling volume and HasOnlySubstanceUnits """
        print("\n### (9) Demo of StochPy's direct method for handling volume and HasOnlySubstanceUnits ###")
        input("press any button to continue\n")
        print(">>> smod.Model('dsmts-02-11.xml.psc')")
        smod.Model('dsmts-001-11.xml.psc')  
        print(">>> smod.DoStochSim(method = 'direct',trajectories=1000,end=50,mode='time')")
        smod.DoStochSim(method = 'direct',trajectories=1000,end=50,mode ='time')
        print(">>> smod.GetRegularGrid()")
        smod.GetRegularGrid()        
        print(">>> smod.PlotAverageSpeciesTimeSeries()")        
        smod.PlotAverageSpeciesTimeSeries()

    def Demo10(self):
        """ Demo of StochPy for doing (preprogrammed) Sequential simulations """
        print("\n### (10) Demo of StochPy for doing (preprogrammed) Sequential simulations ###")
        input("press any button to continue\n")
        print(">>> cmod = stochpy.CellDivision() # start cell division module")
        cmod = stochpy.CellDivision()
        print(">>> cmod.DoCellDivisionStochSim(mode='generations',end = 3, trajectories=1)")
        cmod.DoCellDivisionStochSim(mode='generations',end = 3,trajectories=1)
        print(">>> cmod.PlotSpeciesTimeSeries() # plot time series of species")
        cmod.PlotSpeciesTimeSeries()
        
    def Demo11(self):
        """ StochPy Demo which illustrates delayed SSAs """
        print(">>> smod.Model('Isomerization.psc')")
        smod.Model('Isomerization.psc')
        print(">>> smod.DoStochSim(mode='time',end=10,trajectories=1000)")
        smod.DoStochSim(mode='time',end=10,trajectories=1000)
        print(">>> smod.GetRegularGrid(n_samples=51)")
        smod.GetRegularGrid(n_samples=51)
        print(">>> smod.PlotAverageSpeciesTimeSeries()")
        smod.PlotAverageSpeciesTimeSeries()   
        
        print(">>> smod.SetDelayParameters({'R1':('fixed',5)})")
        smod.SetDelayParameters({'R1':('fixed',5)})
        print(">>> smod.DoDelayedStochSim(mode='time',end=10,trajectories=1000)")
        smod.DoDelayedStochSim(mode='time',end=10,trajectories=1000)
        print(">>> smod.GetRegularGrid()")
        smod.GetRegularGrid()
        print(">>> smod.PlotAverageSpeciesTimeSeries()")
        smod.PlotAverageSpeciesTimeSeries()
    
    def DoDemoSimulations(self):
        self.Demo1()
        self.Demo2()
        self.Demo3()
        self.Demo4()
        self.Demo5()
        self.Demo6()
        self.Demo7()
        self.Demo8()
        self.Demo9()
        self.Demo10()
        self.Demo11()
