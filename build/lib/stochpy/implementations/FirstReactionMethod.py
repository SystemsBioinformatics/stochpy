#! /usr/bin/env python
"""
First Reaction Method
=====================
This module performs the first reaction method Stochastic Simulation Algorithm from Gillespie (1977).

This algorithm is used to generate exact realizations of the Markov jump process. Of course, the algorithm is stochastic, so these realizations are different for each run.

Only molecule populations are specified. Positions and velocities, such as in Molecular Dynamics (MD) are ignored. This makes the algorithm much faster, because non-reactive molecular collisions can be ignored.different
Still, this exact SSA is quite slow, because it insists on simulating every individual reaction event, which takes a lot of time if the reactant population is large.
Furthermore, even larger problems arise if the model contains distinct processes operating on different time-scales.

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: June 23, 2015
"""

from __future__ import division, print_function, absolute_import

############################## IMPORTS ###################################

import sys,copy,time,os
from .StochPyTools import __species__,StochPySSA_Shared,np

########################### END IMPORTS ##################################

class FirstReactionMethod(StochPySSA_Shared):
    """ 
    First Reaction Method from Gillespie (1977)

    This algorithm is used to generate exact realizations of the Markov jump process. Of course, the algorithm is stochastic, so these realizations are   different for each run.
    
    Input:
     - *model_file* filename.psc
     - *model_dir* /home/user/Stochpy/pscmodels/filename.psc
     """
    def __init__(self,model_file,model_dir,IsQuiet=False): 
        self.Parse(model_file,model_dir,IsQuiet=IsQuiet)  
            

    def Execute(self,settings,IsStatusBar=False):       
        """
        Generates T trajectories of the Markov jump process.

        Input:
         - *settings* (class object)   
        """
        if settings.IsSeed:
            np.random.seed(5) 
        
        self._IsInitial = True
        self.settings = settings
        self.sim_t = copy.copy(settings.starttime)
        self.fixed_species_amount = copy.deepcopy(self.parse.fixed_species_amount)
        self.X_matrix = copy.deepcopy(settings.X_matrix)               
        
        #New Volume addition (also modification in self.propensities)
        try:
            self.volume_code = settings.volume_code
        except AttributeError: 
            self.volume_code = "self._current_volume = 1"    # No volume_code present in settings (normal mode)        
          
        self.Propensities()
        
        if not self.sim_t:      
            self.timestep = 1   
            self.sim_output = []
            self.propensities_output = []              
            self.V_output = []        
            self._IsTrackPropensities = copy.copy(settings.IsTrackPropensities)      
            self.SpeciesSelection() 
            self.RateSelection()
            self.SetEvents()     
            if not settings.IsOnlyLastTimepoint:
                self.Initial_Conditions()                            

        nstep_counter = 1
        t1 = time.time()
        while self.sim_t < settings.endtime and self.timestep < settings.timesteps:               
            if self.sim_a_0 <= 0:                         # All reactants got exhausted                 
                 settings.endtime = 10**50
                 break
                 
            #print(self.sim_t)             
            self.RunExactTimestep()                       # Run SSA time step             
            self.HandleEvents()  
            
            # Update Propensities
            if self.sim_t < settings.endtime:
                if not self._IsPerformEvent:
                    self.species_to_update = self.parse.reaction_affects[self.reaction_index] # Determine vars to update                
                else:
                    self.species_to_update  = [s for s in range(self.n_species)]
               
                self.Propensities()
        
                if not settings.IsOnlyLastTimepoint:      # Store Output
                    self.GenerateOutput()
                    
            self._IsPerformEvent = False                  # set to false (or just to make sure).            
            t2 = time.time()        
            if IsStatusBar and t2-t1> 1:
                t1 = time.time()
                sys.stdout.write('\rsimulating {0:s}\r'.format('.'*nstep_counter) ) 
                sys.stdout.flush() 
                nstep_counter+=1
                if nstep_counter > 10:
                    nstep_counter = 1                
                    sys.stdout.write('\rsimulating {0:s}         '.format('.'*nstep_counter))
                    sys.stdout.flush()
        if settings.IsOnlyLastTimepoint or settings.endtime != 10**50:         
           self.GenerateOutput()        
        if IsStatusBar and t1 and not settings.quiet:  
            t1 = time.time()
            sys.stdout.write('\rsimulation done!               \n')   
            

    def RunExactTimestep(self):
        """ Perform a direct SSA time step and pre-generate M random numbers """    
        randoms = np.random.random(self.n_reactions)             # Regenerate for each time step M random numbers
        self.randoms_log = np.log(randoms)*-1     
        #print(self.sim_t,self.sim_a_mu,self.X_matrix)
        self.sim_taus_list = (self.randoms_log[0:self.n_reactions]/self.sim_a_mu).tolist()
        self.sim_tau = min(self.sim_taus_list)                   # Select minimum tau
        
        if (self.sim_t + self.sim_tau) < self.settings.endtime:
            self.sim_t += self.sim_tau            
            self.reaction_index = self.sim_taus_list.index(self.sim_tau)
            try:
                self.X_matrix += self.N_matrix_transpose[self.reaction_index]
                self.timestep += 1
            except MemoryError as ex:
                print(ex)
                sys.exit()            
        else:          
            self.sim_t = self.settings.endtime
            self.reaction_index = np.nan
