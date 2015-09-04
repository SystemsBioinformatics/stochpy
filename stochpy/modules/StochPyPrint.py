 #! /usr/bin/env python
"""
StochPy printing class
======================

Written by T.R. Maarleveld and M. Moinat, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: August 05, 2015
"""

############################ IMPORTS ################################
from __future__ import division, print_function, absolute_import

class PrintingFunctions():

    def PrintSpeciesTimeSeries(self):
        """ Print time simulation output for each generated trajectory """      
        assert self._IsSimulationDone, "First do a stochastic simulation"      
        for n in range(1,self.sim_trajectories_done+1):
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            print('Time', "\t", end="")
            for s_id in self.sim_species_tracked:
                print(s_id,"\t",end="")
            print()
            for timepoint in self.data_stochsim.getSpecies():
                for value in timepoint:
                    print(value,"\t",end="")
                print()            


    def PrintPropensitiesTimeSeries(self):
        """ Print a time series of the propensities each generated trajectory """      
        assert (self._IsTrackPropensities and self._IsSimulationDone), "First do a stochastic simulation with tracking propensities (use the IsTrackPropensities flag in DoStochSim)"
        for n in range(1,self.sim_trajectories_done+1):   	
            if self.sim_trajectories_done > 1:        
                self.GetTrajectoryData(n)
            print('Time', "\t", end="")
            for r_id in self.sim_rates_tracked:
                print(r_id,"\t",end="")
            print()
            for timepoint in self.data_stochsim.getPropensities():
                for value in timepoint:               
                     print(value,"\t",end="")
                print()


    def PrintSpeciesDistributions(self):
        """ Print obtained distributions for each generated trajectory """      
        assert self._IsSimulationDone, "First do a stochastic simulation"
        assert not self._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"
        for n in range(1,self.sim_trajectories_done+1):
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            for i,L_species_dist in enumerate(self.data_stochsim.species_distributions):
                print("Copy number ({0:s})\tPMF".format(self.sim_species_tracked[i]) )
                for m in range(len(L_species_dist[0])):
                    x = L_species_dist[0][m]
                    p_x = L_species_dist[1][m]
                    if not p_x < 0.001:
                        print("{0:d}\t{1:0.3f}".format(x,p_x))
                    else:
                        print("{0:d}\t{1:0.3e}".format(x,p_x))          


    def PrintPropensitiesDistributions(self):
        """ Print obtained distributions for each generated trajectory """     
        assert (self._IsTrackPropensities and self._IsSimulationDone),"First do a stochastic simulation"
        assert not self._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"
        for n in range(1,self.sim_trajectories_done+1):    
            if self.sim_trajectories_done > 1:  
                self.GetTrajectoryData(n)
            
            for j,L_prop_dist in enumerate(self.data_stochsim.propensities_distributions):
                print("Propensity ({0:s})\tPMF".format(self.sim_rates_tracked[j]) )
                for m in range(len(L_prop_dist[0])):                    
                    x = L_prop_dist[0][m]
                    p_x = L_prop_dist[1][m]
                    if not p_x < 0.001:
                        print("{0}\t{1:0.3f}".format(x,p_x))   
                    else:
                        print("{0}\t{1:0.3e}".format(x,p_x))                                           


    def PrintWaitingtimesDistributions(self):
        """ Print obtained waiting times """
        assert not self._IsTauleaping, "Tau-Leaping method does not allow for calculation of waiting times"
        if (not self.data_stochsim.HAS_WAITINGTIMES) and (not self._IsTauleaping):
            self.GetWaitingtimes()
        for n in range(1,self.sim_trajectories_done+1):   
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)           
            for r_id in self.data_stochsim.waiting_times:
                print("Waiting times\t({0:s})".format(r_id) )
                waiting_times_r = self.data_stochsim.waiting_times[r_id]
                for wtime in waiting_times_r:
                    if not wtime < 0.001:
                        print("{0:0.3f}".format(wtime))   
                    else:
                        print("{0:0.3e}".format(wtime))                        
                        

    def PrintSpeciesMeans(self):
        """ Print the means (3 decimals) of each species for the selected trajectory"""
        assert self._IsSimulationDone, "First do a stochastic simulation"
        assert not self._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"
        print("Species\tMean")
        for s_id in self.data_stochsim.species_labels:   
            mu = self.data_stochsim.species_means[s_id]         
            if not mu < 0.001:
                print("{0:s}\t{1:0.3f}".format(s_id,mu))   
            else:
                print("{0:s}\t{1:0.3e}".format(s_id,mu))           


    def PrintSpeciesStandardDeviations(self):
        """ Print the standard deviations (3 decimals) of each species for the selected trajectory"""  
        assert self._IsSimulationDone, "First do a stochastic simulation"
        assert not self._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"           
        print("Species\tStandard Deviation")
        for s_id in self.data_stochsim.species_labels:          
            sigma = self.data_stochsim.species_standard_deviations[s_id]
            if not sigma < 0.001:
                print("{0:s}\t{1:0.3f}".format(s_id,sigma))     
            else:
                print("{0:s}\t{1:0.3e}".format(s_id,sigma))    
                

    def PrintPropensitiesMeans(self): 
        """ Print the means of each propensity for the selected trajectory"""      
        assert (self._IsTrackPropensities and self._IsSimulationDone), "First do a stochastic simulation with tracking propensities (use the IsTrackPropensities flag in DoStochSim)"      
        assert not self._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"
        print("Reaction\tMean")
        for r_id in self.sim_rates_tracked:
            mu = self.data_stochsim.propensities_means[r_id]
            if not mu < 0.001:
                print("{0:s}\t{1:0.3f}".format(r_id,mu))
            else:
                print("{0:s}\t{1:0.3e}".format(r_id,mu))


    def PrintPropensitiesStandardDeviations(self):
        """ Print the standard deviations of each propensity for the selected trajectory"""  
        assert (self._IsTrackPropensities and self._IsSimulationDone), "First do a stochastic simulation with tracking propensities (use the IsTrackPropensities flag in DoStochSim)"    
        assert not self._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"           
        print("Reaction\tStandard Deviation")
        for r_id in self.sim_rates_tracked:            
            std = self.data_stochsim.propensities_standard_deviations[r_id]
            if not std < 0.001:
                print("{0:s}\t{1:0.3f}".format(r_id,std))     
            else:
                print("{0:s}\t{1:0.3e}".format(r_id,std))    


    def PrintWaitingtimesMeans(self):
        """ Print the waiting time means for the selected trajectory """      
        assert not self._IsTauleaping, "Tau-Leaping method does not allow for calculation of waiting times"        
        if (not self.data_stochsim.HAS_WAITINGTIMES) and (not self._IsTauleaping): 
            self.GetWaitingtimes()

        for n in range(1,self.sim_trajectories_done+1):     
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            print("Reaction\tMean")            
            for j,r_id in enumerate(self.SSA.rate_names):                  
                mu = self.data_stochsim.waiting_times_means[j]
                if not mu < 0.001:
                    print("{0:s}\t{1:0.3f}".format(r_id,mu))
                else:
                    print("{0:s}\t{1:0.3e}".format(r_id,mu))              
   

    def PrintWaitingtimesStandardDeviations(self):
        """ Print the waiting time standard deviations for the selected trajectory """
        assert not self._IsTauleaping, "Tau-Leaping method does not allow for calculation of waiting times"    
        if (not self.data_stochsim.HAS_WAITINGTIMES) and (not self._IsTauleaping): 
            self.GetWaitingtimes()
        for n in range(1,self.sim_trajectories_done+1):     
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            print("Reaction\tStandard deviation")            
            for j,r_id in enumerate(self.SSA.rate_names):                             
                std = self.data_stochsim.waiting_times_standard_deviations[j]
                if not std < 0.001:
                    print("{0:s}\t{1:0.3f}".format(s_id,std))     
                else:
                    print("{0:s}\t{1:0.3e}".format(s_id,std))                        


    def PrintAverageSpeciesTimeSeries(self):    
        """ Analyze the average output over all generated trajectories """
        if not self.HAS_AVERAGE:
            print("*** WARNING ***: No regular grid is created yet. Use GetRegularGrid(n_samples) if averaged results are unsatisfactory (e.g. more or less 'points')")
            self.GetRegularGrid()
        for s_id in self.data_stochsim.species_labels:              
            print("\t{0:s} (Mean)\t{0:s} (STD)".format(s_id),end="")
        print()
        for x,t in enumerate(self.data_stochsim_grid.time): 
            print(t,end="")            
            for i in range(len(self.data_stochsim_grid.species_labels)):                                 
                mu = self.data_stochsim_grid.species_means[x,i]
                sigma = self.data_stochsim_grid.species_standard_deviations[x,i]
                if not mu < 0.001 and not sigma < 0.001:
                    print("\t{0:0.3f}\t{1:0.3f}".format(mu,sigma),end="") 
                elif not mu < 0.001:
                    print("\t{0:0.3f}\t{1:0.3e}".format(mu,sigma),end="")
                else:
                    print("\t{0:0.3e}\t{1:0.3e}".format(mu,sigma),end="")                                                     
            print()            
                

    def PrintAveragePropensitiesTimeSeries(self):
        """ Analyze the average output over all generated trajectories """
        assert self._IsTrackPropensities, "First do a stochastic simulation with tracking propensities (use the IsTrackPropensities flag in DoStochSim)"
        if (not self.HAS_AVERAGE) and (self._IsTrackPropensities):
            print("*** WARNING ***: No regular grid is created yet. Use GetRegularGrid(n_samples) if averaged results are unsatisfactory (e.g. more or less 'points')")
            self.GetRegularGrid()
        print("Time",end="")
        for r_id in self.sim_rates_tracked:
            print("\t{0:s} (Mean)\t{0:s} (STD)".format(r_id),end="")
        print()
        for x,t in enumerate(self.data_stochsim_grid.time): 
            print(t,end="")
            for j in range(len(self.sim_rates_tracked)):                 
                mu = self.data_stochsim_grid.propensities_means[x,j]
                sigma = self.data_stochsim_grid.propensities_standard_deviations[x,j]
                if not mu < 0.001 and not sigma < 0.001:
                    print("\t{0:0.3f}\t{1:0.3f}".format(mu,sigma),end="")   
                elif not mu < 0.001:
                    print("\t{0:0.3f}\t{1:0.3e}".format(mu,sigma),end="")   
                else:
                    print("\t{0:0.3e}\t{1:0.3e}".format(mu,sigma),end="")    
            print()
