 #! /usr/bin/env python
"""
StochPyCellDivision plotting class
==================================

Written by T.R. Maarleveld and M. Moinat, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: August 05, 2015
"""
from __future__ import division, print_function, absolute_import

from . import Analysis

import copy

try: import numpy as np
except:
    print("Make sure that the NumPy module is installed")
    print("This program does not work without NumPy")
    print("See http://numpy.scipy.org/ for more information about NumPy")
    sys.exit()
        

class PlottingFunctions():
        
    def PlotSpeciesOverview(self,species2plot = True,sample='extant',main_colors = ['#0072B2','#009E73','#D55E00']): # TODO
        """
        """     
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"   
        assert not self.StochSim._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        assert self._IsAnalyzedExtant , "First analyze the extant cell population"           
        
        if sample not in ['extant','mother','baby']:
            sample = 'extant'
           
        # species selection
        if species2plot == True:
            species_index = 0
            species2plot = self.data_stochsim.species_labels[species_index]
            if not len(self.data_stochsim.species_labels) == 1:
                print("Warning: no species was selected, so the first species '{0}' is selected by default".format(species2plot))
        else:  
            assert species2plot in self.data_stochsim.species_labels,"Species '{0}' is not in the model or species selection".format(species2plot)
            species_index = self.data_stochsim.species_labels.index(species2plot)

        Analysis.plt.figure(self.StochSim.plot.plotnum)    

        # make a grid for plotting
        gs = Analysis.gridspec.GridSpec(1, 3,width_ratios=[0.333,1,0.33333],height_ratios=[1])
        ax1 = Analysis.plt.subplot(gs[2])  

        sample_ = copy.copy(sample)
        if sample.lower() not in ['mother','extant']:
           sample_ ='extant'

        self.PlotSpeciesDistributions(sample=sample_,species2plot=species2plot,title=r'$n_{%s}$' % sample_,xlabel='',IsLegend=False,orientation='horizontal',colors='#E69D00',multiplotting=True)

        max_x = max(self._interdivision_times) # should be fixed, so min should be max
        Analysis.plt.yticks([])

        max_y = np.ceil(max(self.data_stochsim_celldivision.species_at_division[:,species_index]))
        Analysis.plt.ylim([0,max_y])
        
        ax2 =  Analysis.plt.subplot(gs[0]) 

        self.StochSim.plot.plotnum -= 1
        self.PlotSpeciesAtDivisionDistributions(sample=sample,species2plot=species2plot,orientation='horizontal',IsLegend=False,colors='#CC79A7',multiplotting=True)

        self.StochSim.plot.plotnum -= 1
        self.PlotSpeciesAtBirthDistributions(sample=sample,species2plot=species2plot,orientation='horizontal',title=r'$n_{a,%s}$' % sample,IsLegend=False,colors='#56B3E9',multiplotting=True)

        Analysis.plt.ylim([0,max_y])
        Analysis.plt.gca().invert_xaxis()

        ax3 = Analysis.plt.subplot(gs[1]) 
        self.StochSim.plot.plotnum -= 1
        
        # central plot - do a short simulation for a neat time series plot       
        steps_cumulative = self.data_stochsim_celldivision.generation_timesteps.cumsum()
        steps_cumulative = np.insert(steps_cumulative,0,0)    
        main_colors = ['#0072B2','#009E73','#D55E00']
        
        n_generations = 3
        if len(self.data_stochsim_celldivision.generation_timesteps) < n_generations:
            n_generations = len(self.data_stochsim_celldivision.generation_timesteps)
        
        for m in range(1,n_generations+1): 
            n_steps = self.data_stochsim_celldivision.generation_timesteps[-m]        
            current_steps = steps_cumulative[-(m+1)]        
            data = np.column_stack((self.data_stochsim.time[current_steps:current_steps+n_steps],self.data_stochsim.species[current_steps:current_steps+n_steps,species_index]))
            
            data[:,0] = data[:,0] - data[0,0] # reset to zero      
            data = Analysis.getDataForTimeSimPlot(data)
            
            Analysis.plt.plot(data[:,0],data[:,1],color=main_colors[m-1])    
            
        Analysis.plt.title('$n(t)$')
        Analysis.plt.axvline(np.log(2)/self.sim_volume_growth_rate,ls='dashed',color='black',lw=0.5)        
        Analysis.plt.yticks([])
        Analysis.plt.xlim([0,max_x])
        Analysis.plt.ylim([0,max_y])
        self.StochSim.plot.plotnum+=1 
                 
        
    def PlotVolumeOverview(self,sample='extant',main_colors = ['#0072B2','#009E73','#D55E00']): # TODO: multiple traj?
    
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"
        assert not self.StochSim._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        assert self._IsAnalyzedExtant , "First analyze the extant cell population"
        
        if sample not in ['extant','mother','baby']:
            sample = 'extant'
       
        Analysis.plt.figure(self.StochSim.plot.plotnum)
        gs = Analysis.gridspec.GridSpec(2, 3,width_ratios=[0.333,1,0.33333],height_ratios=[1,0.333])
        ax1 = Analysis.plt.subplot(gs[0,2])  

        # do now a long simulation for accurate distributions
        sample_ = copy.copy(sample)
        if sample.lower() not in ['mother','extant']:
           sample_ ='extant'
        self.PlotVolumeDistribution(sample=sample_,orientation='horizontal',colors = '#E69D00',xlabel='',title= r'$V_{%s}$' % sample_,multiplotting=True)

        max_x = max(self._interdivision_times)
        Analysis.plt.yticks([])
        max_y = np.ceil(max(self._volume_at_division))
        Analysis.plt.ylim([0,max_y])

        ax2 =  Analysis.plt.subplot(gs[0,0]) 

        self.StochSim.plot.plotnum -= 1 
        self.PlotVolumeAtDivisionDistribution(sample=sample,orientation='horizontal',title='',xlabel='',ylabel='',colors='#56B3E9')

        self.StochSim.plot.plotnum -= 1
        self.PlotVolumeAtBirthDistribution(sample=sample,orientation='horizontal',title=r'$V_{a,%s}$' % sample,xlabel = 'Volume/cell',colors='#CC79A7',multiplotting=True)

        Analysis.plt.ylim([0,max_y])
        Analysis.plt.gca().invert_xaxis()

        ax3 =  Analysis.plt.subplot(gs[1,1]) 
        self.StochSim.plot.plotnum -= 1

        self.PlotInterdivisionTimeDistribution(sample=sample,bins=25,title=r'$f_{%s}(t)$' % sample,colors='#0072B2',multiplotting=True)
        Analysis.plt.axvline(np.log(2)/self.sim_volume_growth_rate,ls='dashed',color='black',lw=0.5)
        Analysis.plt.gca().invert_yaxis()

        Analysis.plt.xlim([0,max_x])

        ax4 = Analysis.plt.subplot(gs[0,1])  
        self.StochSim.plot.plotnum -= 1

        # central plot - do a short simulation for a neat time series plot       
        steps_cumulative = self.data_stochsim_celldivision.generation_timesteps.cumsum()
        steps_cumulative = np.insert(steps_cumulative,0,0)    
        main_colors = ['#0072B2','#009E73','#D55E00']
        
        n_generations = 3
        if len(self.data_stochsim_celldivision.generation_timesteps) < n_generations:
            n_generations = len(self.data_stochsim_celldivision.generation_timesteps)
        
        for m in range(1,n_generations+1): 
            n_steps = self.data_stochsim_celldivision.generation_timesteps[-m]        
            current_steps = steps_cumulative[-(m+1)]        
            data = np.column_stack((self.data_stochsim.time[current_steps:current_steps+n_steps],self.data_stochsim.volume[current_steps:current_steps+n_steps]))
            
            data[:,0] = data[:,0] - data[0,0] # reset to zero      
            data = Analysis.getDataForTimeSimPlot(data)
            
            Analysis.plt.plot(data[:,0],data[:,1],color=main_colors[m-1])    
            
        Analysis.plt.title('$n(t)$')        
        Analysis.plt.axvline(np.log(2)/self.sim_volume_growth_rate,ls='dashed',color='black',lw=0.5)
        Analysis.plt.xlim([0,max_x])
        Analysis.plt.ylim([0,max_y])
        Analysis.plt.yticks([])
        Analysis.plt.tight_layout()
        self.StochSim.plot.plotnum+=1     
        
        
    def PlotSpeciesDistributions(self, sample='extant',species2plot = True, histtype='step',linestyle = 'solid', linewidth = 1,colors=None,title = 'Species distributions (sample = )',xlabel='Copy number/cell',ylabel='PMF',IsLegend=True,legend_location='upper right',bin_size=1,orientation='vertical',multiplotting=False):
        """
        Plots the species PMF for a sample of extant or mother cells.
        
        Input:
         - *sample* [default = 'extant'] (string) ['extant','mother']
         - *species2plot* [default = True) (list)     
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar]              
         - *linestyle* [default = 'dotted'] (string)         
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = 'Extant Species Distributions'] (string)
         - *xlabel* [default = 'Volume'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *bin_size* [default=None] (integer)
         - *orientation* [default = 'vertical'] (string)
         - *multiplotting* [default = False] (boolean)
        """ 
        assert Analysis._IsPlotting, "Install Matplotlib or use Export2file()"        
            
        species2plot = self.StochSim._getSpecies2Plot(species2plot)        
        species2plot_indices = [self.StochSim.sim_species_tracked.index(d) for d in species2plot]        
        if colors == None:
            colors = self.StochSim.plot.colors

        Analysis.plt.figure(self.StochSim.plot.plotnum)                
        for n in range(1,self.StochSim.sim_trajectories_done+1):
            if self.StochSim.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            
            if sample.lower() == 'mother':    # mother
                distribution2plot = self.data_stochsim.species_distributions
            else:                             # extant
                sample = 'extant'
                assert self._IsAnalyzedExtant, "First analyze the extant cells (.AnalyzeExtantCells())."
                distribution2plot = self.data_stochsim_celldivision.species_extant_distributions        
   
            if '(sample = )' in title:
                title = title.replace('= ','= {0:s}'.format(sample.lower()))
   
            self.StochSim.plot.Distributions(distributions=distribution2plot,datatype= species2plot,labels=self.data_stochsim_celldivision.species_labels,trajectory_index = n-1,linestyle= linestyle,linewidth = linewidth,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,is_legend=IsLegend,legend_location = legend_location,bin_size = bin_size,orientation=orientation,histtype=histtype,multiplotting=multiplotting)  
        self.StochSim.plot.plotnum += 1   
        
    # TODO: multiple trajectories FIX THIS MESH
    def PlotSpeciesAtBirthDistributions(self, sample='extant', bin_size=1, species2plot = True, histtype = 'step', linestyle = 'solid', filled = False, linewidth = 1,colors=None,title = 'Species at birth (sample = )',xlabel='Copy number/cell',ylabel='PMF',IsLegend=True,legend_location='upper right',orientation='vertical',multiplotting=False):
        """ 
        Plots the species at birth distribution for a sample of mother cells 
        
        Input:
         - *sample* [default = 'extant'] (string) ['mother','baby','extant']
         - *bin_size* [default = 1] (integer)
         - *species2plot* [default = True) (list) 
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar] 
         - *linestyle* [default = 'dotted'] (string)
         - *filled* [default = False] (boolean)
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = 'Species at birth (sample = )'] (string)
         - *xlabel* [default = 'Volume'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *orientation* [default = 'vertical'] (string)
         - *multiplotting* [default = False] (boolean)            
        """
        assert len(self._species_at_division) > 0, "Not enough generations to plot species distribution at birth"
        species2plot = self.StochSim._getSpecies2Plot(species2plot)
        
        if len(np.unique(self._interdivision_times)) == 1 and sample.lower() == 'extant': # alle samples are identical, extant is not calculated
            sample = 'baby'
            
        if '(sample = )' in title:
            title = title.replace('= ','= {0:s}'.format(sample.lower()))                       
        
        if sample.lower() in ['mother','baby']: 
            self._PlotSpecies(plottype = "birth",sample=sample,bin_size=bin_size,species2plot=species2plot,histtype=histtype,linestyle=linestyle,filled=filled,linewidth = linewidth,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,orientation=orientation,multiplotting=multiplotting)        
                          
            max_copy_number = 0
            for s_id in species2plot:
                i = self.data_stochsim_celldivision.species_labels.index(s_id)
                max_i = max(self.data_stochsim_celldivision.species_at_division[:,i])
                if max_i > max_copy_number:
                   max_copy_number = max_i        
            if orientation =='vertical':
                Analysis.plt.xlim(xmax=max_i)    
            else:
                Analysis.plt.ylim(ymax=max_i)  
            
        else: # extant,          
            species2plot_indices = [self.StochSim.sim_species_tracked.index(d) for d in species2plot]        
            if colors == None:
                colors = self.StochSim.plot.colors    
            
            Analysis.plt.figure(self.StochSim.plot.plotnum)                
            for n in range(1,self.StochSim.sim_trajectories_done+1):
                if self.StochSim.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                        
                species_distribution = self.data_stochsim_celldivision.species_extant_at_birth_distributions
       
                self.StochSim.plot.Distributions(distributions=species_distribution,datatype= species2plot,labels=self.data_stochsim_celldivision.species_labels,trajectory_index = n-1,linestyle= linestyle,linewidth = linewidth,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,is_legend=IsLegend,legend_location = legend_location,bin_size = bin_size,orientation=orientation,histtype=histtype,multiplotting=multiplotting)  
            self.StochSim.plot.plotnum += 1 
        

    def PlotSpeciesAtDivisionDistributions(self, sample ='extant',bin_size=1, species2plot = True, histtype = 'step', linestyle = 'solid', filled = False, linewidth = 1,colors=None,title = 'Species at division (sample = )',xlabel='Copy number/cell',ylabel='PMF',IsLegend=True,legend_location='upper right',orientation='vertical',multiplotting=False):
        """ 
        Plots the species at division distribution for a sample of mother cells
        
        Input:      
         - *sample* [default = 'extant'] (string) ['extant','mother']
         - *bin_size* [default = 1] (integer)
         - *species2plot* [default = True) (list) 
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar] 
         - *linestyle* [default = 'dotted'] (string)
         - *filled* [default = False] (boolean)
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = 'Volume'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *orientation* [default = 'vertical'] (string)  
         - *multiplotting* [default = False] (boolean)                   
        """        
        assert len(self._species_at_division) > 0, "Not enough generations to plot a mother species distribution"        
        if len(np.unique(self._interdivision_times)) == 1 and sample.lower() != 'mother': # all samples are identical, extant is not calculated
            sample = 'mother'        
            
        if '(sample = )' in title:
            title = title.replace('= ','= {0:s}'.format(sample.lower()))            
        
        if sample.lower() == 'mother':            
            self._PlotSpecies(plottype = "division",sample='mother',bin_size=bin_size,species2plot=species2plot,histtype=histtype,linestyle=linestyle,filled=filled, linewidth = linewidth,colors=colors,title =title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location, orientation=orientation, multiplotting=multiplotting)
        else:  # baby, or extant
            species2plot = self.StochSim._getSpecies2Plot(species2plot)        
            species2plot_indices = [self.StochSim.sim_species_tracked.index(d) for d in species2plot]        
            if colors == None:
                colors = self.StochSim.plot.colors    
            
            Analysis.plt.figure(self.StochSim.plot.plotnum)                
            for n in range(1,self.StochSim.sim_trajectories_done+1):
                if self.StochSim.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                        
                if sample.lower() == 'baby':
                    species_distribution = self.data_stochsim_celldivision.species_baby_at_division_distributions
                else:
                    species_distribution = self.data_stochsim_celldivision.species_extant_at_division_distributions
       
                self.StochSim.plot.Distributions(distributions=species_distribution,datatype= species2plot,labels=self.data_stochsim_celldivision.species_labels,trajectory_index = n-1,linestyle= linestyle,linewidth = linewidth,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,is_legend=IsLegend,legend_location = legend_location,bin_size = bin_size,orientation=orientation,histtype=histtype,multiplotting=multiplotting)  
            self.StochSim.plot.plotnum += 1                               
             

    def PlotVolumeDistribution(self, sample='extant', histtype='step',linestyle = 'solid', linewidth = 1,colors=None,title = 'Volume distribution (sample = )', xlabel='Volume',ylabel='PDF',IsLegend=True,legend_location='upper right',orientation='vertical',multiplotting=False):
        """
        Plots the volume PMF for a sample of extant or mother cells.
        
        We use the number of bins that were used to generate the volume distribution
        
        Input:      
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar]              
         - *linestyle* [default = 'dotted'] (string)         
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = 'Extant Cell Volume'] (string)
         - *xlabel* [default = 'Volume'] (string)
         - *ylabel* [default = 'PDF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)       
         - *orientation* [default = 'vertical'] (string)
         - *multiplotting* [default = False] (boolean)         
        """ 
        assert Analysis._IsPlotting, "Install Matplotlib or use Export2file()"
       
        if colors == None:
            colors =  ['#0000FF','#00CC00','#FF0033','#FF00CC','#6600FF','#FFFF00','#000000','#CCCCCC',
                       '#00CCFF','#99CC33','#FF6666','#FF99CC','#CC6600','#003300','#CCFFFF','#9900FF','#CC6633','#FFD700','#C0C0C0']   
        elif isinstance(colors,str):
            colors = [colors]    
                    
        j=0
        Analysis.plt.figure(self.StochSim.plot.plotnum)                        
        for n in range(1,self.StochSim.sim_trajectories_done+1):
            if self.StochSim.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)

            if j >= len(colors): # make sure that plotting works also if not enough colors are provided
                j=0

            if sample.lower() == 'mother':                         
                output = Analysis.plt.hist(self.data_stochsim.volume,bins=20,color = colors[j],histtype=histtype,normed=True,orientation=orientation) # TODO hardcoded # bins
            else:
                assert self._IsAnalyzedExtant, "First analyze the extant cells (.AnalyzeExtantCells())."          
                sample = 'extant'
                x = self.data_stochsim_celldivision.volume_extant_distribution[0]
                y = self.data_stochsim_celldivision.volume_extant_distribution[1]                
                output = Analysis.plt.hist(x,len(x)-1,weights=y,normed=True,ls=linestyle,lw=linewidth,histtype=histtype,color=colors[j],orientation=orientation,align='left')
            j+=1

        if '(sample = )' in title:
            title = title.replace('= ','= {0:s}'.format(sample.lower()))

        Analysis.plt.title(title)
        if orientation.lower() == 'horizontal':
            Analysis.plt.xlabel(ylabel)
            Analysis.plt.ylabel(xlabel)    
            if multiplotting:
                Analysis.plt.xticks([0,max(output[0])*1.2])
        else:             
            Analysis.plt.xlabel(xlabel)
            Analysis.plt.ylabel(ylabel)    
            if multiplotting:
                Analysis.plt.yticks([0,max(output[0])*1.2])
        self.StochSim.plot.plotnum+=1         

   
    def PlotVolumeAtDivisionDistribution(self,sample='extant', bins=10, histtype = 'step', linestyle = 'solid', filled = False, linewidth = 1,colors=None,title = 'Volume distribution at division (sample = )',xlabel='Volume',ylabel='PDF',orientation='vertical',multiplotting=False):
        """
        Plot the volume distribution of mother cells
        
        Input:
         - *sample* [default = 'extant'] (string) ['mother','extant']        
         - *bin_size* [default = 0.1] (float)
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar] 
         - *linestyle* [default = 'dotted'] (string)
         - *filled* [default = False] (boolean)
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = 'Volume'] (string)
         - *ylabel* [default = 'PDF'] (string)
         - *orientation* [default = 'vertical'] (string)
         - *multiplotting* [default = False] (boolean)
        """
        assert len(self._volume_at_division) > 0, "Not enough generations to plot a mother volume distribution"
        if len(np.unique(self._interdivision_times)) == 1 and sample.lower() != 'mother': # all samples are identical, extant is not calculated
            sample = 'mother'      

        self._PlotVolume(plottype = "division",sample=sample,bins= bins, histtype = histtype, linestyle = linestyle, filled = filled, linewidth = linewidth,colors=colors,title = title,xlabel=xlabel,ylabel=ylabel,orientation=orientation,multiplotting=multiplotting)


    def PlotVolumeAtBirthDistribution(self,sample='extant',bins= 10, histtype = 'step', linestyle = 'solid', filled = False, linewidth = 1,colors=None,title = 'Volume distribution at birth (sample = )',xlabel='Volume',ylabel='PDF',orientation='vertical',multiplotting=False):
        """
        Plot the volume distribution of daughter cells
        
        Input:
         - *sample* [default = 'extant'] (string) ['mother','baby','extant']
         - *bins* [default = 0.1] (float)
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar] 
         - *linestyle* [default = 'dotted'] (string)
         - *filled* [default = False] (boolean)
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = 'Volume'] (string)
         - *ylabel* [default = 'PDF'] (string)
         - *orientation* [default = 'vertical'] (string)
         - *multiplotting* [default = False] (boolean)         
        """
        assert len(self._volume_at_birth) > 0, "Not enough generations to plot a daughter volume distribution"        
        
        if len(np.unique(self._interdivision_times)) == 1 and sample.lower() == 'extant': # all samples are identical, extant is not calculated
            sample = 'mother'
        
        self._PlotVolume(plottype = "birth",sample = sample,bins= bins, histtype = histtype, linestyle = linestyle, filled = filled, linewidth = linewidth,colors=colors,title = title,xlabel=xlabel,ylabel=ylabel,orientation=orientation,multiplotting=multiplotting)        


    def PlotInterdivisionTimeDistribution(self,sample='extant', bins=10, histtype = 'step', linestyle = 'solid', filled = False, linewidth = 1,colors=None,title = False,xlabel='Interdivision Time',ylabel='PDF',IsLegend=True,legend_location='upper right',multiplotting=False):
        """
        Plots Interdivision Time distribution (between cell divisions) for a certain *trajectories* or for 'all' trajectories together        
        
        Input:
         - *sample* [default = 'extant'] (string) ['extant','baby','mother']
         - *bins* [default = 10] (integer)
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar] 
         - *linestyle* [default = 'dotted'] (string)
         - *filled* [default = False] (boolean)
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = False] (string)
         - *xlabel* [default = 'Interdivision Time'] (string)
         - *ylabel* [default = 'PDF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *multiplotting* [default=False] (boolean)
        """         
        assert len(self._interdivision_times) > 1, "Not enough generations to plot an interdivision time distribution"        
        assert Analysis._IsPlotting, "Install Matplotlib or use Export2file()"
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"
        
        if colors == None:
            colors =  ['#0000FF','#00CC00','#FF0033','#FF00CC','#6600FF','#FFFF00','#000000','#CCCCCC',
                       '#00CCFF','#99CC33','#FF6666','#FF99CC','#CC6600','#003300','#CCFFFF','#9900FF','#CC6633','#FFD700','#C0C0C0']           
        elif isinstance(colors,str):
            colors = [colors]    
        
        Analysis.plt.figure(self.StochSim.plot.plotnum)    
        j=0    
        for n in range(1,self.StochSim.sim_trajectories_done+1):   
            if self.StochSim.sim_trajectories > 1:
                self.GetTrajectoryData(n)
            
            f_m, tau = np.histogram(self.data_stochsim_celldivision.interdivision_times, bins=bins,density=True)
            tau = np.array([(x+y)/2. for x,y in zip(tau,tau[1:])]) # center stuff

            f_b = 0.5 * np.exp(self.sim_volume_growth_rate*tau)*f_m
            f_e = f_b * 2*(1-np.exp(-self.sim_volume_growth_rate*tau))

            if j >= len(colors):                  
                j=0               

            if sample.lower() == 'baby':             
                Analysis.plt.plot(tau,f_b,drawstyle='steps-mid',color=colors[j]) 
                max_y = max(f_b)               
                if title == False:
                    title = r'$f_b(t)$'
            elif sample.lower() == 'mother':
                Analysis.plt.plot(tau,f_m,drawstyle='steps-mid',color=colors[j])
                max_y = max(f_m)  
                if title == False:
                    title = r'$f_m(t)$'
            else: # extant
                Analysis.plt.plot(tau,f_e,drawstyle='steps-mid',color=colors[j])
                max_y = max(f_e)  
                if title == False:
                    title = r'$f_e(t)$'                
            j+=1                    
 
        Analysis.plt.xlim(xmin=0)
        Analysis.plt.ylim(ymin=0)
        
        self.StochSim.plot.plotnum += 1
        
        Analysis.plt.title(title)
        Analysis.plt.xlabel(xlabel)
        Analysis.plt.ylabel(ylabel)        
       
        if multiplotting:
            Analysis.plt.yticks([0,round(max_y*1.2,2)])     
             
        
    def PlotCellAgeDistribution(self, bins=50,  histtype = 'step', linestyle = 'solid', filled = False, linewidth = 1,colors=None,title = '',xlabel='Cell Age',ylabel='PDF',orientation='vertical'):
        """
        Plot the cell age distribution
        
        Input:
         - *bins [default = 50] (integer)
         - *histtype* [default = 'step'] (string) ['step','stepfilled','bar] 
         - *linestyle* [default = 'dotted'] (string)
         - *filled* [default = False] (boolean)
         - *linewidth* [default = 1] (float)
         - *colors* [default = None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = 'Volume'] (string)
         - *ylabel* [default = 'PDF'] (string)
         - *orientation* [default = 'vertical'] (string)
        """        
        # TODO: for different samples?
        self._PlotVolume(plottype = "cellage",sample='mother',bins= bins, histtype = histtype, linestyle = linestyle, filled = filled, linewidth = linewidth,colors=colors,title = title,xlabel=xlabel,ylabel=ylabel,orientation=orientation)
        

    def _PlotSpecies(self, plottype, sample='mother', bin_size=1, species2plot = True, histtype = 'step', linestyle = 'solid', filled = False, linewidth = 1,colors=None,title = '',xlabel='Copy number/cell',ylabel='PMF',IsLegend=True,legend_location='upper right',orientation='vertical',multiplotting=False):
        """ 
        Normed=False because the total probability is determined by summation not by integration.
        
        *** For internal use only ***
        """
        assert Analysis._IsPlotting, "Install Matplotlib or use Export2file()"
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"
            
        species2plot = self.StochSim._getSpecies2Plot(species2plot)
        
        species2plot_indices = [self.StochSim.sim_species_tracked.index(d) for d in species2plot]        
        if colors == None:
            colors =  ['#0000FF','#00CC00','#FF0033','#FF00CC','#6600FF','#FFFF00','#000000','#CCCCCC',
                       '#00CCFF','#99CC33','#FF6666','#FF99CC','#CC6600','#003300','#CCFFFF','#9900FF','#CC6633','#FFD700','#C0C0C0']           
        elif isinstance(colors,str):
            colors = [colors]              
                
        Analysis.plt.figure(self.StochSim.plot.plotnum)
        for n in range(1,self.StochSim.sim_trajectories_done+1):   
            if self.StochSim.sim_trajectories > 1:
                self.GetTrajectoryData(n)
            
            if plottype.lower() == 'division':
                data = self.data_stochsim_celldivision.species_at_division
            elif plottype.lower() == 'birth':   
                if sample.lower() == 'mother':
                    data = self.data_stochsim_celldivision.species_at_birth # mother cells
                elif sample.lower() =='baby':                 
                    data = np.vstack((self.data_stochsim_celldivision.species_at_birth,self.data_stochsim_celldivision.species_at_birth_not_tracked)) # baby cells
    
            for m,s_index in enumerate(species2plot_indices):
                if len(species2plot_indices) == 1:
                    m = n-1                
                dat_min = data[:,s_index].min() 
                dat_max = data[:,s_index].max()
                n_bins = 1 + (dat_max-dat_min) /bin_size # Just take one trajectory as reference
                L_bin_edges = np.linspace(dat_min-bin_size/2.0,dat_max+bin_size/2.0,n_bins+1) 
                     
                output = Analysis.plt.hist(data[:,s_index], L_bin_edges, normed=True, histtype = histtype, ls = linestyle, fill = filled, lw = linewidth, color = colors[m],orientation=orientation) 
        
        self.StochSim.plot.plotnum += 1        
        Analysis.plt.title(title)
        if orientation.lower() == 'horizontal':
            Analysis.plt.xlabel(ylabel)
            Analysis.plt.ylabel(xlabel) 
            Analysis.plt.xlim(xmin=0)
            Analysis.plt.ylim(ymin=0) 
            if multiplotting:
                Analysis.plt.xticks([0,max(output[0])*1.2])                      
        else:             
            Analysis.plt.xlabel(xlabel)
            Analysis.plt.ylabel(ylabel)   
            Analysis.plt.xlim(xmin=0)
            Analysis.plt.ylim(ymin=0)       
            if multiplotting:
                Analysis.plt.yticks([0,max(output[0])*1.2])               
        if IsLegend:
            Analysis.plt.legend(species2plot,numpoints=1,frameon=True,loc=legend_location)        
        
            
    def _PlotVolume(self, plottype, sample = 'extant',bins=10, histtype = 'step', linestyle = 'solid', filled = False, linewidth = 1,colors=None,title = '',xlabel='Volume',ylabel='PDF',orientation='vertical',multiplotting=False):
        """ 
        Here, we can directly use the hist() function. This is not possible for e.g. time series data, where we have to take into account the irregular time points between consecutive firings of reactions. Species distributions can therefore not be directly be generated with the hist() function.
        
        continuous in time: cell age distribution, interdivision time, mother cell volume, daughter cell volume. 
        Normed=True because the total probability is determined by integration.
        
        *** For internal use only ***
        """
        assert Analysis._IsPlotting, "Install Matplotlib or use Export2file()"
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"
        
        if colors == None:
            colors =  ['#0000FF','#00CC00','#FF0033','#FF00CC','#6600FF','#FFFF00','#000000','#CCCCCC',
                       '#00CCFF','#99CC33','#FF6666','#FF99CC','#CC6600','#003300','#CCFFFF','#9900FF','#CC6633','#FFD700','#C0C0C0']           
        elif isinstance(colors,str):
            colors = [colors]
            
        if '(sample = )' in title:
            title = title.replace('= ','= {0:s}'.format(sample.lower()))               
        
        Analysis.plt.figure(self.StochSim.plot.plotnum)    
        j=0    
        for n in range(1,self.StochSim.sim_trajectories_done+1):   
            if self.StochSim.sim_trajectories > 1:
                self.GetTrajectoryData(n)
            if sample.lower() in ['mother','baby']:    
                # We now use mid as alignment for the data
                if plottype.lower() == 'division':
                    data = self.data_stochsim_celldivision.volume_at_division
                    align='mid'
                elif plottype.lower() == 'birth':
                    if sample.lower() == 'baby':    
                        data = self.data_stochsim_celldivision.volume_at_birth_all # sample of baby's
                        align='mid'
                    elif sample.lower() == 'mother':                
                        data = self.data_stochsim_celldivision.volume_at_birth # sample of mothers
                        align='mid'
                elif plottype.lower() in ['age','ages','cellage','cellages']:
                    data = self.data_stochsim_celldivision.ages_deterministic
                    align = 'mid'         
                    
                if j >= len(colors):                  
                    j=0    
                
                output = Analysis.plt.hist(data,bins= bins,normed=True,histtype=histtype,ls=linestyle,fill=filled,lw=linewidth,color=colors[j],align=align,orientation=orientation)
            
            else: # extant
                assert self._IsAnalyzedExtant, "First analyze the extant cells (.AnalyzeExtantCells())."    
                if plottype.lower() == 'division':
                    x = self.data_stochsim_celldivision.volume_extant_at_division_distribution[0]
                    y = self.data_stochsim_celldivision.volume_extant_at_division_distribution[1]                
                elif plottype.lower() == 'birth':
                    x = self.data_stochsim_celldivision.volume_extant_at_birth_distribution[0]
                    y = self.data_stochsim_celldivision.volume_extant_at_birth_distribution[1]

                if j >= len(colors):                  
                    j=0   
                
                output = Analysis.plt.hist(x,len(x)-1,weights = y,ls = linestyle,lw = linewidth,histtype = histtype,color = colors[j],orientation=orientation,align='left',normed=True)
            
            j+=1
        
        Analysis.plt.xlim(xmin=0)
        Analysis.plt.ylim(ymin=0)   
        
        self.StochSim.plot.plotnum += 1
        Analysis.plt.title(title)
        if orientation.lower() == 'horizontal':
            Analysis.plt.xlabel(ylabel)
            Analysis.plt.ylabel(xlabel)    
            if multiplotting:
                Analysis.plt.xticks([0,max(output[0])*1.2])
        else:             
            Analysis.plt.xlabel(xlabel)
            Analysis.plt.ylabel(ylabel)    
            if multiplotting:
                Analysis.plt.yticks([0,max(output[0])*1.2])
                

    def PlotVolumeTimeSeries(self,n_events2plot = 10000,linestyle = 'solid', linewidth = 1,marker = '',colors = None,title = '', xlabel = 'Time', ylabel = 'Volume',deterministic=True):
        """
        StochPy Volume Time Series Plot
        
        Input:
         - *n_events2plot* [default = 10000] (integer)       
         - *linestyle* [default = 'solid'] (string)       
         - *linewidth* [default = 1] (float)
         - *marker* [default = ''] (string)
         - *colors* [default = None] (list)
         - *title* [default = ''] (string)
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Volume'] (string)   
         - *deterministic* [default = True] (boolean)     
        """
        assert Analysis._IsPlotting, "Install matplotlib or use Export2file()"
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"        
        if str(n_events2plot).lower() == 'all':
            n_events2plot = self.data_stochsim.simulation_timesteps       
                
        Analysis.plt.figure(self.StochSim.plot.plotnum)        
        if colors == None: 
            colors = self.StochSim.plot.colors    
        elif isinstance(colors,str):
            colors = [colors]    
        j=0
        for n in range(1,self.StochSim.sim_trajectories_done+1):   
            if self.StochSim.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
                
            if j >= len(colors): # make sure that plotting works also if not enough colors are provided
                j=0
                                
            if deterministic:                
                Arr_volume = Analysis.getDataForTimeSimPlot(self.data_stochsim_celldivision.getVolumeDeterministic(),quiet=self.StochSim._IsQuiet)
            else:
                Arr_volume = Analysis.getDataForTimeSimPlot(self.data_stochsim.getVolume(),quiet= self.StochSim._IsQuiet)
            
            Analysis.plt.plot(Arr_volume[:,0],Arr_volume[:,1], marker, ls = linestyle, lw = linewidth, color = colors[j])  
            j+=1 
        Analysis.plt.title(title)
        Analysis.plt.xlabel(xlabel)
        Analysis.plt.ylabel(ylabel)
        self.StochSim.plot.plotnum+=1     
               

    def PlotSpeciesTimeSeries(self,plottype = 'copy numbers',n_events2plot = 10000,species2plot = True,linestyle = 'solid',linewidth = 1,marker = '',colors = None,title = '',xlabel='Time',ylabel='Copy number/cell',IsLegend=True,legend_location='upper right'):
        """
        Plot time simulation output for each generated trajectory.

        Input:
         - *plottype* [default = 'copy numbers'] (string) ['concentrations','copy numbers']
         - *n_events2plot* [default = 10000] (integer)
         - *species2plot* [default = True] as a list ['S1','S2'] 
         - *linestyle* [default = 'solid'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = ''] ('v','o','s',',','*','.')
         - *colors* [default = None] (list)
         - *title* [default = '']  (string)
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Copy number/cell'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string)
        """     
        assert Analysis._IsPlotting, "Install matplotlib or use Export2file()"
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"
        assert not self.StochSim._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        species2plot = copy.copy(self.StochSim._getSpecies2Plot(species2plot))
        if str(n_events2plot).lower() == 'all':
            n_events2plot = self.data_stochsim.simulation_timesteps
                
        for n in range(1,self.StochSim.sim_trajectories_done+1):   
            if self.StochSim.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            
            if plottype.lower() in ['concentration','concentrations']:               
                plot_data = self.data_stochsim.getSpeciesConcentrations()
                species_labels = copy.copy(self.data_stochsim.species_concentrations_labels)
                for s in self._species_extracellular:
                    if s in species2plot:
                         species2plot.remove(s)                                
                if ylabel == 'Copy number/cell': 
                     ylabel = 'Concentration'           
            elif plottype.lower() in ['amounts','amount','copy numbers','copy_number']: 
                plot_data = self.data_stochsim.getSpecies()                
                species_labels = copy.copy(self.data_stochsim.species_labels)
            else:
                raise Warning("{0} is not a valid plottype argument. Valid plottypes are 'concentrations' and 'copy numbers'".format(plottype) )            
       
            self.StochSim.plot.TimeSeries(plot_data,n_events2plot,species2plot,species_labels,n-1,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location)           
        self.StochSim.plot.plotnum += 1
        
    
    def PlotSpeciesVolumeTimeSeries(self,plottype ='copy numbers', n_events2plot = 10000,species2plot = True,linestyle = 'solid',linewidth = 1, marker = '',colors = None,title= 'Species and Volume Time Series',xlabel = 'Time',ylabel = 'Copy number/cell',IsLegend = True,legend_location=1):
        """
        Plot time simulation output for each generated trajectory        

        Input:
         - *plottype* [default = 'copy numbers'] (string) ['concentrations','copy numbers']
         - *n_events2plot* [default = 10000] (integer)
         - *species2plot* [default = True] as a list ['S1','S2'] 
         - *linestyle* [default = 'solid'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = ''] ('v','o','s',',','*','.')
         - *colors* [default = None] (list)
         - *title* [default = 'Species and Volume Time Series']  (string)    
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Copy number/cell'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string)   
        """
        assert Analysis._IsPlotting, "Install matplotlib or use Export2file()"
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"
        assert not self.StochSim._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        #Make two subplots
        gs = Analysis.gridspec.GridSpec(2,1)
        plotnum = self.StochSim.plot.plotnum
        Analysis.plt.figure(plotnum);
        
        #Plot species in upper panel
        ax1 = Analysis.plt.subplot(gs[0])
        self.StochSim.plot.plotnum = plotnum
        self.PlotSpeciesTimeSeries(plottype,n_events2plot, species2plot, linestyle, linewidth, marker, colors,title,xlabel,ylabel,IsLegend=IsLegend,legend_location=legend_location)
        Analysis.plt.xlabel('')
        
        #Plot volume in lower panel
        ax2 = Analysis.plt.subplot(gs[1])
        self.StochSim.plot.plotnum = plotnum
        self.PlotVolumeTimeSeries(n_events2plot,linestyle, linewidth, marker, colors,title='',xlabel=xlabel)
        
                
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
         - *legend_location* [default = 'upper right'] (string)   
        """
        self.StochSim.PlotPropensitiesTimeSeries(n_events2plot=n_events2plot,rates2plot=rates2plot,linestyle=linestyle,linewidth=linewidth,marker=marker,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location)
        
        
    def PlotPropensitiesDistributions(self,rates2plot = True, linestyle = 'solid',linewidth = 1,colors=None,title ='', xlabel='Propensity',ylabel='PMF',IsLegend=True,legend_location = 'upper right',bin_size=1):
        """
        Plots the probability mass function for each generated trajectory

        Input:
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] (string)
         - *linewidth* [default = 1] (float)
         - *colors* (list)
         - *title* [default = ''] (string)     
         - *xlabel* [default = 'Propensity'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *bin_size* [default=1] (integer)
        """
        self.StochSim.PlotPropensitiesDistributions(rates2plot=rates2plot,linestyle=linestyle,linewidth=linewidth,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,bin_size=bin_size)  
         

    def PlotWaitingtimesDistributions(self,rates2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',colors = None,title = '',xlabel=r'inter-event time $t$',ylabel='PDF',IsLegend=True,legend_location='upper right'):
        """
        Plot event waiting time distributions
        
        default: PlotWaitingtimesDistributions() plots waiting times for all rates
      
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
         - *legend_location* [default = 'upper right'] (string)   
        """     
        self.StochSim.PlotWaitingtimesDistributions(rates2plot=rates2plot,linestyle=linestyle,linewidth=linewidth,marker=marker,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location)


    def PlotAverageVolumeTimeSeries(self,linestyle = '', linewidth = 1,marker = 'o',color = 'blue',title = 'Average Volume Time Series (# of trajectories = )', xlabel = 'Time', ylabel = 'Volume',nstd=1):
        """
        StochPy Average Volume Time Series Plot
        
        Input:             
         - *linestyle* [default = 'solid'] (string)       
         - *linewidth* [default = 1] (float)
         - *marker* [default = ''] (string)
         - *color* (string)
         - *title* [default = 'Average Volume Time Series (# of trajectories = )'] (string)
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Volume'] (string)   
         - *nstd* [default = 1] (float)     
        """
        assert Analysis._IsPlotting, "Install matplotlib or use Export2file()"
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"        
        if not self.StochSim.HAS_AVERAGE: 
            print("*** WARNING ***: No regular grid is created yet. Use GetRegularGrid(n_samples) if averaged results are unsatisfactory (e.g. more or less 'samples')")
            self.GetRegularGrid()  
        
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.sim_trajectories_done))    
        
        plotnum = self.StochSim.plot.plotnum
        Analysis.plt.figure(plotnum);
            
        Analysis.plt.errorbar(self.data_stochsim_grid.time,self.data_stochsim_grid.volume_means[0],yerr = nstd*self.data_stochsim_grid.volume_standard_deviations[0],ls=linestyle,lw=linewidth,marker=marker,color=color)

        Analysis.plt.title(title)
        Analysis.plt.xlabel(xlabel)
        Analysis.plt.ylabel(ylabel)
        self.StochSim.plot.plotnum+=1 
        
        
    def PlotAverageSpeciesVolumeTimeSeries(self,plottype ='copy numbers', species2plot = True,linestyle = 'solid',linewidth = 1, marker = '',colors = None,title= 'Average Species and Volume Time Series (# of trajectories = )',xlabel = 'Time', ylabel = 'Copy number/cell', IsLegend = True,legend_location=1,nstd=1): 
        """
        Plot time simulation output for each generated trajectory        

        Input:
         - *plottype* [default = 'copy numbers'] (string) ['concentrations','copy numbers']
         - *species2plot* [default = True] as a list ['S1','S2'] 
         - *linestyle* [default = 'solid'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = ''] ('v','o','s',',','*','.')
         - *colors* [default = None] (list)
         - *title* [default = 'Species Time Series and Volume Plot']  (string)   
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Copy number/cell'] (string)             
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string) 
         - *nstd* [default = 1] (float)              
        """
        assert Analysis._IsPlotting, "Install matplotlib or use Export2file()"
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"   
        assert not self.StochSim._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        #Make two subplots
        gs = Analysis.gridspec.GridSpec(2,1)
        plotnum = self.StochSim.plot.plotnum
        Analysis.plt.figure(plotnum);
        
        #Plot species in upper panel
        ax1 = Analysis.plt.subplot(gs[0])
        self.StochSim.plot.plotnum = plotnum
        self.PlotAverageSpeciesTimeSeries(plottype,species2plot, linestyle, linewidth, marker, colors,title,xlabel,ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)
        Analysis.plt.xlabel('')
        
        #Plot volume in lower panel
        ax2 = Analysis.plt.subplot(gs[1])
        self.StochSim.plot.plotnum = plotnum
        self.PlotAverageVolumeTimeSeries(linestyle, linewidth, marker, title='',nstd=nstd)   
        
        
    def PlotAverageSpeciesExtantDistributions(self,species2plot = True,linestyle = 'None',linewidth = 1,marker = 'o',colors = None,title = 'Average Extant Species Distributions (# of trajectories = )',xlabel='Copy number/cell',ylabel='PMF',IsLegend=True,legend_location = 'upper right',nstd=1): 
        """
        Plot the average species distributions For each species Amount, the mean and standard deviation are plotted      

        Input:
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Extant Species Distributions (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Copy number/cell'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *nstd* [default=1] (float)
        """
        assert Analysis._IsPlotting, "Install Matplotlib or use Export2file()"
        assert self._IsAnalyzedExtant, "First analyze the simulations (.AnalyzeExtantCells())."      
        if not self.data_stochsim_grid.HAS_AVERAGE_SPECIES_EXTANT_DISTRIBUTIONS:
            self.GetAverageSpeciesExtantDistributions()
        
        species2plot = self.StochSim._getSpecies2Plot(species2plot)                
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.StochSim.sim_trajectories_done))

        self.StochSim.plot.AverageDistributions(self.data_stochsim_grid.species_extant_distributions_means,self.data_stochsim_grid.species_extant_distributions_standard_deviations,nstd,species2plot,self.StochSim.sim_species_tracked,linestyle,linewidth,marker,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.StochSim.plot.plotnum+=1        
                   
        
    def PlotAverageSpeciesTimeSeries(self,plottype='copy numbers',species2plot = True,linestyle = 'None',linewidth = 1,marker = 'o',markersize=5,colors = None,title = 'Average Species Time Series (# of trajectories = )',xlabel='Time',ylabel='Copy number/cell',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted 
        
        Input:
         - *plottype* [default = 'copy numbers'] (string) ['concentrations','copy numbers']
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *markersize* [default = 5] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Time Series (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Time'] (string)
         - *ylabel* [default = 'Copy number/cell'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string)   
         - *nstd* [default=1] (float)
        """ 
        assert Analysis._IsPlotting, "Install matplotlib or use Export2file()"  
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"   
        assert not self.StochSim._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        
        if not self.StochSim.HAS_AVERAGE: 
            print("*** WARNING ***: No regular grid is created yet. Use GetRegularGrid(n_samples) if averaged results are unsatisfactory (e.g. more or less 'samples')")
            self.GetRegularGrid()                      
        species2plot = copy.copy(self.StochSim._getSpecies2Plot(species2plot))        
        if '(# of trajectories = )' in title:
            title = title.replace('= ','= {0:d}'.format(self.StochSim.sim_trajectories_done))            
            
        if plottype.lower() in ['concentration','concentrations']:               
            means = copy.copy(self.data_stochsim_grid.species_concentrations_means)
            stds = copy.copy(self.data_stochsim_grid.species_concentrations_standard_deviations)
            for s in self._species_extracellular:
                if s in species2plot:
                     species2plot.remove(s)                                
            if ylabel == 'Copy number/cell': 
                 ylabel = 'Concentration'           
        elif plottype.lower() in ['amounts','amount','copy numbers','copy_number']: 
            means = copy.copy(self.data_stochsim_grid.species_means)
            stds = copy.copy(self.data_stochsim_grid.species_standard_deviations)        
        
        self.StochSim.plot.AverageTimeSeries(means,stds,self.data_stochsim_grid.time,nstd,species2plot,self.StochSim.sim_species_tracked,linestyle,linewidth,marker,markersize,colors,title,xlabel,ylabel,IsLegend,legend_location)
        self.StochSim.plot.plotnum+=1       
        
        
    def PlotAveragePropensitiesTimeSeries(self,rates2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',markersize=5,colors = None,title = 'Average Propensities Time Series (# of trajectories = )',xlabel='Time',ylabel='Propensity',IsLegend=True,legend_location='upper right',nstd=1): 
        """        
        Plot the average propensities For each time point, the mean and standard deviation are plotted 
        
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
         - *legend_location* [default = 'upper right'] (string)   
         - *nstd* [default=1] (float)
        """              
        self.StochSim.PlotAveragePropensitiesTimeSeries(rates2plot=rates2plot,linestyle=linestyle,linewidth=linewidth,marker=marker,markersize=markersize,colors=colors,title=title,xlabel=xlabel,ylabel=ylabel,legend_location=legend_location,nstd=nstd)          
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 
        
        
    def PlotAverageSpeciesDistributions(self,species2plot = True,linestyle = 'None',linewidth = 1,marker = 'o',colors = None,title = 'Average Species Distributions (# of trajectories = )',xlabel='Copy number/cell',ylabel='PMF',IsLegend=True,legend_location = 'upper right',nstd=1): 
        """
        Plot the average species distributions For each species Amount, the mean and standard deviation are plotted      

        Input:
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Distributions (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Copy number/cell'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string)   
         - *nstd* [default=1] (float)
        """    
        self.StochSim.PlotAverageSpeciesDistributions(species2plot = species2plot,linestyle = linestyle,linewidth = linewidth,marker = marker,colors = colors,title = title, xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 
        
        
    def PlotAverageSpeciesDistributionsConfidenceIntervals(self,species2plot=True,colors = None,title = 'Average Species Distributions (# of trajectories = )',xlabel='Copy number/cell',ylabel='PMF',IsLegend=True,legend_location = 'upper right',nstd=1):
        """
        Plot the average species distributions For each species Amount, the mean and standard deviation are plotted      

        Input:
         - *species2plot* [default = True] as a list ['S1','S2']
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Distributions (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Copy number/cell'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)      
         - *legend_location* [default = 'upper right'] (string)   
         - *nstd* [default=1] (float)
        """
        self.StochSim.PlotAverageSpeciesDistributionsConfidenceIntervals(species2plot=species2plot,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 
        

    def PlotAveragePropensitiesDistributions(self,rates2plot = True,linestyle = 'None',linewidth = 1,marker = 'o',colors = None,title = 'Average Propensities Distributions (# of trajectories = )',xlabel='Propensity',ylabel='PMF',IsLegend=True,legend_location = 'upper right',nstd=1): 
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
         - *legend_location* [default = 'upper right'] (string)   
         - *nstd* [default=1] (float)
        """  
        self.StochSim.PlotAveragePropensitiesDistributions(rates2plot = rates2plot,linestyle = linestyle,linewidth = linewidth,marker = marker,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid)
         
        
    def PlotAveragePropensitiesDistributionsConfidenceIntervals(self,rates2plot = True,colors = None,title = 'Average Propensities Distributions (# of trajectories = )',xlabel='Propensity',ylabel='PMF',IsLegend=True,legend_location='upper right',nstd=1):        
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted      

        Input:
         - *rates2plot* [default = True] as a list ['R1','R2']       
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Propensities Distributions (# of trajectories = ... )' ] (string)
         - *xlabel* [default = 'Propensity'] (string)
         - *ylabel* [default = 'PMF'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string)   
         - *nstd* [default=1] (float)
        """
        self.StochSim.PlotAveragePropensitiesDistributionsConfidenceIntervals(rates2plot = rates2plot,colors = colors,title = title, xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid)        
        
        
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
        self.StochSim.PlotSpeciesAutocorrelations(nlags = nlags,species2plot=species2plot,linestyle = linestyle,linewidth = linewidth, marker = marker,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 


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
        self.StochSim.PlotSpeciesAutocovariances(nlags = nlags,species2plot=species2plot,linestyle = linestyle,linewidth = linewidth, marker = marker,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 
        
        
    def PlotAverageSpeciesAutocorrelations(self,nlags=-1,species2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',markersize=5,colors = None,title = 'Average Species Autocorrelations (# of trajectories = )',xlabel=r'Lag ($\tau$)',ylabel='Auto-correlation',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted       

        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)      
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *markersize* [default = 5] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Autocorrelations (# of trajectories = ... )' ] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-correlation'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """        
        self.StochSim.PlotAverageSpeciesAutocorrelations(nlags = nlags,species2plot=species2plot,linestyle = linestyle,linewidth = linewidth, marker = marker,markersize=markersize,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 
        

    def PlotAverageSpeciesAutocovariances(self,nlags=-1,species2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',markersize=5,colors = None,title = 'Average Species Autocovariances (# of trajectories = )',xlabel=r'Lag ($\tau$)',ylabel='Auto-covariance',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average time simulation result. For each time point, the mean and standard deviation are plotted       

        Input:

         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *species2plot* [default = True] as a list ['S1','S2']
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = 'o'] ('v','o','s',',','*','.')
         - *markersize* [default = 5] (float)
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Species Autocovariances (# of trajectories = ... )' ] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-covariance'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        self.StochSim.PlotAverageSpeciesAutocovariances(nlags = nlags,species2plot=species2plot,linestyle = linestyle,linewidth = linewidth, marker = marker,markersize=markersize,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 
        
        
    def PlotPropensitiesAutocorrelations(self,nlags=-1,rates2plot=True,linestyle = 'None', linewidth = 1, marker = 'o',colors = None,title = '',xlabel=r'Lag ($\tau$)',ylabel='Auto-correlation',IsLegend=True,legend_location='upper right'):
        """
        PlotPropensitiesAutocorrelations(rates2plot=True,linestyle = 'None',linewidth = 1, marker = 'o',colors = None,title = '',xlabel=r'Lag ($\tau$)',ylabel='Auto-correlation')
        
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
        self.StochSim.PlotPropensitiesAutocorrelations(nlags = nlags,rates2plot=rates2plot,linestyle = linestyle,linewidth = linewidth, marker = marker,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid)         
        
        
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
        self.StochSim.PlotPropensitiesAutocovariances(nlags = nlags,rates2plot=rates2plot,linestyle = linestyle,linewidth = linewidth, marker = marker,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 
        
        
    def PlotAveragePropensitiesAutocorrelations(self,nlags=-1,rates2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',colors = None,title = 'Average Propensities Autocorrelations (# of trajectories = )',xlabel=r'Lag ($\tau$)',ylabel='Auto-correlation',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average propensities autocorrelation result for different lags. For each lag, the mean and standard deviation are plotted       

        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *rates2plot* [default = True] as a list ['R1','R2']      
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *marker* [default = ','] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default = 'Average Propensities Autocorrelations (# of trajectories = ... )'] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-correlation'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """  
        self.StochSim.PlotAveragePropensitiesAutocorrelations(nlags = nlags,rates2plot=rates2plot,linestyle = linestyle,linewidth = linewidth, marker = marker,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)      
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid)   
        
        
    def PlotAveragePropensitiesAutocovariances(self,nlags=-1,rates2plot = True,linestyle = 'None',linewidth = 1, marker = 'o',colors = None,title = 'Propensities Autocovariances (# of trajectories = )',xlabel=r'Lag ($\tau$)',ylabel='Auto-covariance',IsLegend=True,legend_location='upper right',nstd=1): 
        """
        Plot the average propensities autocorrelation result for different lags. For each lag, the mean and standard deviation are plotted

        Input:
         - *nlags* [default = -1] (integer) 1,2,3 ... -1 where 3 means calculate the autocorrelation for the first 3 lags and -1 for all lags  
         - *rates2plot* [default = True] as a list ['R1','R2']      
         - *linestyle* [default = 'dotted'] dashed, solid, and dash_dot (string)
         - *linewidth* [default = 1] (float)
         - *marker* [default = ','] ('v','o','s',',','*','.')
         - *colors* [default =  None] (list)
         - *title* [default =  'Average Propensities Autocovariances (# of trajectories = ... ) ] (string)
         - *xlabel* [default = r'Lag ($\tau$)'] (string)
         - *ylabel* [default = 'Auto-covariance'] (string)
         - *IsLegend* [default = True] (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """          
        self.StochSim.PlotAveragePropensitiesAutocovariances(nlags = nlags,rates2plot=rates2plot,linestyle = linestyle,linewidth = linewidth, marker = marker,colors = colors,title = title,xlabel=xlabel,ylabel=ylabel,IsLegend=IsLegend,legend_location=legend_location,nstd=nstd)        
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid) 
