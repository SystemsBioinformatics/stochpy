 #! /usr/bin/env python
"""
StochPy plotting class
======================

Written by T.R. Maarleveld and M. Moinat, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: August 05, 2015
"""

from __future__ import division, print_function, absolute_import

import stochpy

from . import Analysis

try: 
    import numpy as np  
except ImportError:
    print("Make sure that the NumPy module is installed")
    print("This program does not work without NumPy")
    print("See http://numpy.scipy.org/ for more information about NumPy")
    sys.exit()

class PlottingFunctions():

    def PlotSpeciesTimeSeries(self,n_events2plot = 10000,species2plot = True,linestyle = 'solid',linewidth = 1,marker = '',colors = None,title = '',xlabel='Time',ylabel='Copy Number',IsLegend=True,legend_location='upper right'):
        """
        Plot time simulation output for each generated trajectory

        Input:
         - *n_events2plot* [default = 10000] (integer)
         - *species2plot* [default = True] as a list ['S1','S2'] 
         - *linestyle* [default = 'solid'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = ''] ('v','o','s',',','*','.')
         - *colors* [default = None] (list)
         - *title* [default = '']  (string)
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Copy Number'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """    
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"
        assert self._IsSimulationDone, "First do a stochastic simulation"
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
           
        species2plot = self._getSpecies2Plot(species2plot)
                      
        if str(n_events2plot).lower() == 'all':
                n_events2plot = self.data_stochsim.simulation_timesteps        
        for n in range(1,self.sim_trajectories_done+1):    
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            self.plot.TimeSeries(self.data_stochsim.getSpecies(),n_events2plot,species2plot,self.data_stochsim.species_labels,n-1,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location) # Plot time sim
        self.plot.plotnum+=1        


    def PlotPropensitiesTimeSeries(self,n_events2plot = 10000,rates2plot = True,linestyle = 'solid',linewidth = 1,marker = '',colors = None,title = '',xlabel='Time',ylabel='Propensity',IsLegend=True,legend_location='upper right'):
        """
        Plot time simulation output for each generated trajectory

        Input:
         - *n_events2plot* [default = 10000] (integer)
         - *rates2plot* [default = True]: species as a list ['S1','S2']
         - *marker* [default = ''] ('v','o','s',',','*','.')
         - *linestyle* [default = 'solid']: dashed, dotted, and solid (string)
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Propensity'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"
        assert self._IsSimulationDone and self._IsTrackPropensities, "First do a stochastic simulation with tracking propensities (use the IsTrackPropensities flag in DoStochSim)"
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"        
        
        rates2plot = self._getRates2Plot(rates2plot)    
        
        if str(n_events2plot).lower() == 'all':
            n_events2plot = self.data_stochsim.simulation_timesteps
            
        for n in range(1,self.sim_trajectories_done+1): 
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            self.plot.TimeSeries(self.data_stochsim.getPropensities(),n_events2plot,rates2plot,self.sim_rates_tracked,n-1,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location)            
        self.plot.plotnum+=1        


    def PlotSpeciesDistributions(self,species2plot = True, histtype = 'step',linestyle = 'solid',linewidth = 1,colors=None,title = '',xlabel='Copy number', ylabel='PMF',IsLegend=True,legend_location='upper right',bin_size=1,orientation='vertical'):  
        """       
        Plots the PMF for each generated trajectory    

        Input:
         - *species2plot* [default = True] as a list ['S1','S2']
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar]
         - *linestyle* [default = 'dotted'] (string)
         - *linewidth* [default = 1] (float)
         - *colors* (list)
         - *title* [default = ''] (string)     
         - *xlabel* [default = 'Copy number'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *bin_size* [default=None] (integer)        
         - *orientation* [default = 'vertical'] (string)
        """       
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"
        assert self._IsSimulationDone, "First do a stochastic simulation"
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        species2plot = self._getSpecies2Plot(species2plot)  
        
        for n in range(1,self.sim_trajectories_done+1):   
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            self.plot.Distributions(distributions =self.data_stochsim.species_distributions,datatype=species2plot,labels=self.data_stochsim.species_labels,trajectory_index=n-1,linestyle=linestyle,linewidth=linewidth,
            colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,is_legend=IsLegend,legend_location = legend_location,bin_size= bin_size,histtype = histtype,orientation = orientation)                    
        self.plot.plotnum += 1        
                

    def PlotPropensitiesDistributions(self,rates2plot = True,histtype='step', linestyle = 'solid',linewidth = 1,colors=None,title='',xlabel='Propensity', ylabel='PMF', IsLegend=True, legend_location='upper right',bin_size=1,orientation = 'vertical'):
        """
        Plots the PMF for each generated trajectory    
  
        Input:
         - *rates2plot* [default = True] as a list ['R1','R2']
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar]
         - *linestyle* [default = 'dotted'] (string)
         - *linewidth* [default = 1] (float)
         - *colors* (list)
         - *title* [default = ''] (string)     
         - *xlabel* [default = 'Propensity'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *bin_size* [default=1] (integer)
         - *orientation* [default = 'vertical'] (string)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"
        assert self._IsSimulationDone and self._IsTrackPropensities, "First do a stochastic simulation with tracking propensities (use the IsTrackPropensities flag in DoStochSim)"  
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        rates2plot = self._getRates2Plot(rates2plot)    
        
        for n in range(1,self.sim_trajectories_done+1):    
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)            
           
            self.plot.Distributions(distributions= self.data_stochsim.propensities_distributions,datatype = rates2plot,labels = self.sim_rates_tracked,trajectory_index =n-1,linestyle=linestyle,linewidth=linewidth,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,is_legend=IsLegend,legend_location=legend_location,bin_size= bin_size,histtype=histtype,orientation=orientation)
        self.plot.plotnum += 1
        
        
            

    def PlotWaitingtimesDistributions(self,rates2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',colors = None,title = '',xlabel=r'inter-event time $t$',ylabel='PDF',IsLegend=True,legend_location='upper right'):
        """    
        Plot event waiting time distributions            

        Input:
         - *rates2plot* [default = True]  as a list of strings ["R1","R2"]
         - *linestyle* [default = 'None'] dashed, dotted, dash_dot, and solid (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = 'inter-event time t'] (string)
         - *ylabel* [default = 'PDF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """    
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"
        assert not self._IsTauleaping, "Tau-Leaping method does not allow for calculation of waiting times"
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        if (not self.data_stochsim.HAS_WAITINGTIMES) and (not self._IsTauleaping):
            self.GetWaitingtimes()
            
        waitingtime_r_ids = list(self.data_stochsim.waiting_times)
        if rates2plot == True:
            rates2plot = np.sort(waitingtime_r_ids)
        else:
            rates2plot = self._getRates2Plot(rates2plot)
            for r_id in rates2plot:
                if r_id + '_Completion' in waitingtime_r_ids: # Automatically also select completion waiting times
                    rates2plot.append( r_id + '_Completion' )                             
        
        for n in range(1,self.sim_trajectories_done+1): 
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            self.plot.WaitingtimesDistributions(self.data_stochsim.waiting_times,rates2plot,n-1,linestyle,linewidth, marker,colors,title,xlabel,ylabel,IsLegend,legend_location)            
        self.plot.plotnum+=1
        
        
    def PlotAverageSpeciesTimeSeries(self,species2plot = True,linestyle = 'None',linewidth = 1,marker = 'o',markersize = 5,colors = None,title = 'Average Species Time Series (# of trajectories = )',xlabel='Time',ylabel='Copy number',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted 
        
        Input:
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *markersize* [default = 5] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Time (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Copy number'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"  
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        if not self.HAS_AVERAGE: 
            print("*** WARNING ***: No regular grid is created yet. Use GetRegularGrid(n_samples) if averaged results are unsatisfactory (e.g. more or less 'samples')")
            self.GetRegularGrid()        
        
        species2plot = self._getSpecies2Plot(species2plot)
                
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))
            
        self.plot.AverageTimeSeries(self.data_stochsim_grid.species_means,self.data_stochsim_grid.species_standard_deviations,self.data_stochsim_grid.time,nstd,species2plot,self.sim_species_tracked,linestyle,linewidth,marker,markersize,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.plot.plotnum+=1
            

    def PlotAverageSpeciesDistributions(self,species2plot = True,linestyle = 'None',linewidth = 1,marker = 'o',colors = None,title = 'Average Species Distributions (# of trajectories = )',xlabel='Copy number',ylabel='PMF',IsLegend=True,legend_location = 'upper right',nstd=1): 
        """
        Plot the average species distributions For each species Amount, the mean and standard deviation are plotted      

        Input:
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Distributions (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Copy number'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"  
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"        
        if not self.data_stochsim_grid.HAS_AVERAGE_SPECIES_DISTRIBUTIONS:
            self.GetAverageSpeciesDistributions()
        
        species2plot = self._getSpecies2Plot(species2plot)
                
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))

        self.plot.AverageDistributions(self.data_stochsim_grid.species_distributions_means,self.data_stochsim_grid.species_distributions_standard_deviations,nstd,species2plot,self.sim_species_tracked,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.plot.plotnum+=1
             

    def PlotAverageSpeciesDistributionsConfidenceIntervals(self,species2plot=True,colors = None,title = 'Average Species Distributions (# of trajectories = )',xlabel='Copy number',ylabel='PMF',IsLegend=True,legend_location='upper right',nstd=1):
        """
        Plot the average species distributions For each species Amount, the mean and standard deviation are plotted      

        Input:
         - *species2plot* [default = True] as a list ['S1','S2']
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Distributions (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Copy number'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)   
         - *legend_location* [default = 'upper right'] (string/integer)   
         - *nstd* [default=1] (float)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"  
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        if not self.data_stochsim_grid.HAS_AVERAGE_SPECIES_DISTRIBUTIONS:
            self.GetAverageSpeciesDistributions()
        
        species2plot = self._getSpecies2Plot(species2plot)
        
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))
            
        self.plot.AverageDistributionsCI(self.data_stochsim_grid.species_distributions_means,self.data_stochsim_grid.species_distributions_standard_deviations,nstd,species2plot,self.sim_species_tracked,colors,title,xlabel,ylabel,IsLegend,legend_location)    
        self.plot.plotnum+=1
             

    def PlotAveragePropensitiesDistributionsConfidenceIntervals(self,rates2plot = True,colors = None,title = 'Average Propensities Distributions (# of trajectories = )',xlabel='Propensity',ylabel='PMF',IsLegend=True,legend_location='upper right',nstd=1):
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted      

        Input:
         - *rates2plot* [default = True] as a list ['R1','R2']       
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Time (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Propensity'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"  
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        if not self.data_stochsim_grid.HAS_AVERAGE_PROPENSITIES_DISTRIBUTIONS:
            self.GetAveragePropensitiesDistributions()          

        rates2plot = self._getRates2Plot(rates2plot)    
        
        self.plot.AverageDistributionsCI(self.data_stochsim_grid.propensities_distributions_means,self.data_stochsim_grid.propensities_distributions_standard_deviations,nstd,rates2plot,self.sim_rates_tracked,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.plot.plotnum+=1
                

    def PlotAveragePropensitiesDistributions(self,rates2plot = True,linestyle = 'None',linewidth = 1,marker = 'o',colors = None,title = 'Average Propensities Distributions (# of trajectories = )',xlabel='Propensity',ylabel='PMF',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted      

        Input:
         - *rates2plot* [default = True] as a list ['R1','R2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Propensities Distributions (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Propensity'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)
        """  
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"  
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        if not self.data_stochsim_grid.HAS_AVERAGE_PROPENSITIES_DISTRIBUTIONS:
            self.GetAveragePropensitiesDistributions()
         
        rates2plot = self._getRates2Plot(rates2plot)    
                
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))

        self.plot.AverageDistributions(self.data_stochsim_grid.propensities_distributions_means,self.data_stochsim_grid.propensities_distributions_standard_deviations,nstd,rates2plot,self.sim_rates_tracked,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.plot.plotnum+=1 
            

    def PlotAveragePropensitiesTimeSeries(self,rates2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',markersize=5,colors = None,title = 'Average Propensities Time Series (# of trajectories = )',xlabel='Time',ylabel='Propensity',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average propensities. For each time point, the mean and standard deviation are plotted 
        
        Input:
         - *rates2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *markersize* [default = 5] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Propensities Time Series (# of trajectories = ...)' ] (string)
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Propensity'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)
        """      
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"       
        assert self._IsTrackPropensities, "First do a stochastic simulation with tracking propensities (use the IsTrackPropensities flag in DoStochSim)"
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        if (not self.HAS_AVERAGE) and (self._IsTrackPropensities): 
            print("*** WARNING ***: No regular grid is created yet. Use GetRegularGrid(n_samples) if averaged results are unsatisfactory)")
            self.GetRegularGrid()      
            
        rates2plot = self._getRates2Plot(rates2plot)  
                
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))

        self.plot.AverageTimeSeries(self.data_stochsim_grid.propensities_means,self.data_stochsim_grid.propensities_standard_deviations,                   self.data_stochsim_grid.time,nstd,rates2plot,self.sim_rates_tracked,linestyle,linewidth,marker,markersize,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.plot.plotnum+=1          
        
    def PlotSpeciesAutocorrelations(self,nlags = -1,species2plot=True,linestyle = 'None',linewidth = 1,marker = 'o',colors = None,title = '',xlabel=r'Lag ($\tau$)',ylabel='Auto-correlation',IsLegend=True,legend_location='upper right'):
        """
        Plot species autocorrelations
        
        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-correlation'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"     
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        if not self.data_stochsim_grid.HAS_SPECIES_AUTOCORRELATIONS:
            print("*** WARNING ***: Autocorrelations are not yet calculated. StochPy automatically calculates autocorrelations with pre-defined settings. You can use GetSpeciesAutocorrelations(species2calc=True,n_samples=51)")
            self.GetSpeciesAutocorrelations(species2calc = species2plot)
        
        species2plot = self._getSpecies2Plot(species2plot)
                
        for n in range(self.sim_trajectories_done):       
            self.plot.Autocorrelations(self.data_stochsim_grid.time[:nlags],self.data_stochsim_grid.species_autocorrelations[:,n],   species2plot,self.sim_species_tracked,n,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location)
            
        self.plot.plotnum+=1
        

    def PlotSpeciesAutocovariances(self,nlags = -1,species2plot=True,linestyle = 'None', linewidth = 1, marker = 'o',colors = None,title = '',xlabel=r'Lag ($\tau$)',ylabel='Auto-covariance',IsLegend=True,legend_location='upper right'):
        """
        Plot species auto-covariances
        
        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-covariance'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"     
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        if not self.data_stochsim_grid.HAS_SPECIES_AUTOCOVARIANCES:
            print("*** WARNING ***: Autocovariances are not yet calculated. StochPy automatically calculates autocovariances with pre-defined settings. You can use GetSpeciesAutocovariances(species2calc=True,n_samples=51)")
            self.GetSpeciesAutocovariances(species2calc = species2plot)
        
        species2plot = self._getSpecies2Plot(species2plot)
            
        for n in range(self.sim_trajectories_done):     
            self.plot.Autocorrelations(self.data_stochsim_grid.time[:nlags],self.data_stochsim_grid.species_autocovariances[:,n],     species2plot,self.sim_species_tracked,n,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location)            
        self.plot.plotnum+=1
           

    def PlotAverageSpeciesAutocorrelations(self,nlags=-1,species2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',markersize=5,colors = None,title = 'Average Species Autocorrelations(# of trajectories = )',xlabel=r'Lag ($\tau$)',ylabel='Auto-correlation',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted       

        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)      
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *markersize* [default = 1] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Autocorrelations(# of trajectories = ... )' ] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-correlation'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)   
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"     
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        if not self.data_stochsim_grid.HAS_SPECIES_AUTOCORRELATIONS:
            print("*** WARNING ***: Autocorrelations are not yet calculated. StochPy automatically calculates autocorrelations with pre-defined settings. You can use GetSpeciesAutocorrelations(species2calc=True,n_samples=51)")
            self.GetSpeciesAutocorrelations(species2calc = species2plot)
      
        species2plot = self._getSpecies2Plot(species2plot)    
            
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))
        
        self.plot.AverageTimeSeries(self.data_stochsim_grid.species_autocorrelations_means[:nlags][:],self.data_stochsim_grid.species_autocorrelations_standard_deviations[:nlags][:],self.data_stochsim_grid.time[:nlags],nstd,species2plot,self.sim_species_tracked, linestyle,linewidth,marker,markersize,colors, title,xlabel,ylabel,IsLegend, legend_location)
        self.plot.plotnum+=1
        

    def PlotAverageSpeciesAutocovariances(self,nlags=-1,species2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',markersize=5,colors = None,title = 'Average Species Autocovariances (# of trajectories = )',xlabel=r'Lag ($\tau$)',ylabel='Auto-covariance',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted       

        Input:

         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *markersize* [default = 1] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Autocovariances (# of trajectories = ... )' ] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-covariance'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
             
        if not self.data_stochsim_grid.HAS_SPECIES_AUTOCOVARIANCES:
            print("*** WARNING ***: Autocovariances are not yet calculated. StochPy automatically calculates autocovariances with pre-defined settings. You can use GetSpeciesAutocovariances(species2calc=True,n_samples=51)")
            self.GetSpeciesAutocovariances(species2calc = species2plot)
      
        species2plot = self._getSpecies2Plot(species2plot)
            
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))
        
        self.plot.AverageTimeSeries(self.data_stochsim_grid.species_autocovariances_means[:nlags][:],self.data_stochsim_grid.species_autocovariances_standard_deviations[:nlags][:],self.data_stochsim_grid.time[:nlags],nstd,species2plot,self.sim_species_tracked,linestyle,linewidth, marker, markersize, colors, title, xlabel, ylabel, IsLegend, legend_location)
        self.plot.plotnum+=1
            

    def PlotPropensitiesAutocorrelations(self,nlags=-1,rates2plot=True,linestyle = 'None', linewidth = 1, marker = 'o',colors = None,title = '',xlabel=r'Lag ($\tau$)',ylabel='Auto-correlation',IsLegend=True,legend_location='upper right'):
        """
        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *rates2plot* [default = True] as a list ['R1','R2']      
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = ','] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-correlation'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"     
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        if not self.data_stochsim_grid.HAS_PROPENSITIES_AUTOCORRELATIONS:
            print("*** WARNING ***: Autocorrelations are not yet calculated. StochPy automatically calculates autocorrelations with pre-defined settings. You can use GetPropensitiesAutocorrelations(rates=True,n_samples=51)")
            self.GetPropensitiesAutocorrelations(rates = rates2plot,n_samples=self.data_stochsim.simulation_endtime)
             
        rates2plot = self._getRates2Plot(rates2plot)    
                    
        for n in range(self.sim_trajectories_done):     
            self.plot.Autocorrelations(self.data_stochsim_grid.time[:nlags],self.data_stochsim_grid.propensities_autocorrelations[:,n],rates2plot,self.sim_rates_tracked,n,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.plot.plotnum+=1  
        

    def PlotPropensitiesAutocovariances(self,nlags=-1,rates2plot=True,linestyle = 'None',linewidth = 1, marker = 'o',colors = None,title = '',xlabel=r'Lag ($\tau$)',ylabel='Auto-covariance',IsLegend=True,legend_location='upper right'):
        """
        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *rates2plot* [default = True] as a list ['R1','R2']      
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = ','] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-covariance'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"     
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        if not self.data_stochsim_grid.HAS_PROPENSITIES_AUTOCOVARIANCES:
            print("*** WARNING ***: Autocovariances are not yet calculated. StochPy automatically calculates autocovariances with pre-defined settings. You can use GetPropensitiesAutocovariances(rates=True,n_samples=51)")
            self.GetPropensitiesAutocovariances(rates = rates2plot,n_samples=self.data_stochsim.simulation_endtime)
             
        rates2plot = self._getRates2Plot(rates2plot)    
                    
        for n in range(self.sim_trajectories_done):     
            self.plot.Autocorrelations(self.data_stochsim_grid.time[:nlags],self.data_stochsim_grid.propensities_autocovariances[:,n],rates2plot,self.sim_rates_tracked,n,linestyle,linewidth, marker,colors,title,xlabel,ylabel,IsLegend,legend_location)                        
        self.plot.plotnum+=1  


    def PlotAveragePropensitiesAutocorrelations(self,nlags=-1,rates2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',markersize=5,colors = None,title = 'Average Propensities Autocorrelations(# of trajectories = )',xlabel=r'Lag ($\tau$)',ylabel='Auto-correlation',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average propensities autocorrelation result for different lags. For each lag, the mean and standard deviation are plotted       

        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *rates2plot* [default = True] as a list ['R1','R2']      
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *marker* [default = ','] ('v','o','s',',','*','.')
         - *markersize* [default = 1] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Propensities Autocorrelation(# of trajectories = ... )' ] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-correlation'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)         
        """  
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"      
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
           
        if not self.data_stochsim_grid.HAS_PROPENSITIES_AUTOCORRELATIONS:
            print("*** WARNING ***: Autocorrelations are not yet calculated. StochPy automatically calculates autocorrelations with pre-defined settings. You can use GetPropensitiesAutocorrelations(rates=True,n_samples=51)")
            self.GetPropensitiesAutocorrelations(rates = rates2plot,n_samples=self.data_stochsim.simulation_endtime)
        
        rates2plot = self._getRates2Plot(rates2plot)
                    
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done)) 

        self.plot.AverageTimeSeries(self.data_stochsim_grid.propensities_autocorrelations_means[:nlags][:], self.data_stochsim_grid.propensities_autocorrelations_standard_deviations[:nlags][:], self.data_stochsim_grid.time[:nlags], nstd, rates2plot, self.sim_rates_tracked, linestyle,linewidth,marker,markersize,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.plot.plotnum+=1


    def PlotAveragePropensitiesAutocovariances(self,nlags=-1,rates2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',markersize=5,colors = None,title = 'Average Propensities Autocovariances (# of trajectories = )',xlabel=r'Lag ($\tau$)',ylabel='Auto-covariance',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average propensities autocorrelation result for different lags. For each lag, the mean and standard deviation are plotted

        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *rates2plot* [default = True] as a list ['R1','R2']      
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = ','] ('v','o','s',',','*','.')
         - *markersize* [default = 1] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Propensities Autocovariances (# of trajectories = )' ] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-covariance'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)         
        """  
        assert stochpy._IsPlotting, "Install matplotlib or use Export2file()"         
        assert not self._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        if not self.data_stochsim_grid.HAS_PROPENSITIES_AUTOCOVARIANCES:
            print("*** WARNING ***: Autocovariances are not yet calculated. StochPy automatically calculates autocovariances with pre-defined settings. You can use GetPropensitiesAutocovariances(rates=True,n_samples=51)")
            self.GetPropensitiesAutocovariances(rates = rates2plot,n_samples=self.data_stochsim.simulation_endtime)
        
        rates2plot = self._getRates2Plot(rates2plot)
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))
        self.plot.AverageTimeSeries(self.data_stochsim_grid.propensities_autocovariances_means[:nlags][:],self.data_stochsim_grid.propensities_autocovariances_standard_deviations[:nlags][:], self.data_stochsim_grid.time[:nlags],nstd,rates2plot, self.sim_rates_tracked, linestyle,linewidth, marker,markersize,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.plot.plotnum+=1               
