#! /usr/bin/env python
"""
Direct Method
=============
This module performs the Direct Stochastic Simulation Algorithm from Gillespie (1977) [1].

This algorithm is used to generate exact realizations of the Markov jump process. Of course, the algorithm is stochastic, so these realizations are different for each run.

Only molecule populations are specified. Positions and velocities, such as in Molecular Dynamics (MD) are ignored. This makes the algorithm much faster, because non-reactive molecular collisions can be ignored.different
Still, this exact SSA is quite slow, because it insists on simulating every individual reaction event, which takes a lot of time if the reactant population is large. Furthermore, even larger problems arise if the model contains distinct processes operating on different time scales [2].

[1] Gillespie D.T (1977), "Exact stochastic simulation of coupled chemical reactions", J.Phys. Chem. 81:2340-2361
[2] Wilkinson D.J (2009), "Stochastic Modelling for quantitative description of heterogeneous biological systems", Nat Rev Genet; 0(2):122-133 

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: June 23, 2015
"""

from __future__ import division, print_function, absolute_import

__doc__ = """
          Direct Method
          =============
          This program performs the direct Stochastic Simulation Algorithm from Gillespie (1977) [1].This algorithm is used to generate exact realizations of the Markov jump process. Of course, the algorithm is stochastic, so these realizations are different for each run.
          
          Only molecule populations are specified. Positions and velocities, such as in Molecular Dynamics (MD) are ignored. This makes the algorithm much faster, because non-reactive molecular collisions can be ignored.
          Still, this exact SSA is quite slow, because it insists on simulating every individual reaction event, which takes a lot of time if the reactant population is large. Furthermore, even larger problems arise if the model contains distinct processes operating on different time-scales [2].
          
          [1] Gillespie D.T (1977), "Exact stochastic simulation of coupled chemical reactions", J.Phys. Chem. 81:2340-2361
          [2] Wilkinson D.J (2009), "Stochastic Modelling for quantitative description of heterogeneous biological systems", Nat Rev Genet; 0(2):122-133
          """
############################# IMPORTS ####################################

import sys,copy,time,os,operator
from .StochPyTools import __species__,StochPySSA_Shared,np
    
########################### END IMPORTS ##################################
  
class DirectMethod(StochPySSA_Shared):
    """ 
    Direct Stochastic Simulation Algorithm from Gillespie (1977) [1].

    This algorithm is used to generate exact realizations of the Markov jump process. Of course, the algorithm is stochastic, so these realizations are different for each run.

    [1] Gillespie D.T (1977), "Exact stochastic simulation of coupled chemical reactions", J.Phys. Chem. 81:2340-2361

    Input:  
     - *model_file* filename.psc
     - *model_dir* /home/user/Stochpy/pscmodels/filename.psc   
    """
    def __init__(self,model_file,model_dir,IsQuiet=False):    
        self.Parse(model_file,model_dir,IsQuiet=IsQuiet)   
             

    def Execute(self,settings,IsStatusBar=False):       
        """
        Generates a trajectory of the Markov jump process.

        Input:
         - *settings* (class object)   
        """  
        if settings.IsSeed:
            np.random.seed(5)     
        
        self._IsInitial = True
        self.settings = settings
        self.sim_t = copy.copy(settings.starttime)  # does not have to start at zero if we perform sequential simulations         
        self.X_matrix = copy.deepcopy(settings.X_matrix) 
        self.fixed_species_amount = copy.deepcopy(self.parse.fixed_species_amount)                  
        
        try:
            self.volume_code = settings.volume_code
        except AttributeError:                      # No volume_code present in settings
            self.volume_code = "self._current_volume = 1"     
        
        #self.species_to_update = [s for s in range(self.n_species)] # ensure that the first run updates all species   
        self.Propensities() 
            
        if not self.sim_t: 
            self.timestep = 1   
            self.sim_output = []
            self.propensities_output = []              
            self.V_output = []
            self._IsTrackPropensities = copy.copy(settings.IsTrackPropensities)
            self.SpeciesSelection()   
            self.RateSelection()                           
            self.SetEvents()  # April 15, moved into here, because otherwise each new cell division cycle starts with a time event, if specified            
            if not settings.IsOnlyLastTimepoint:
                self.Initial_Conditions()                            
                            
        nstep_counter = 1
        t1 = time.time()
        while (self.sim_t < settings.endtime) and (self.timestep < settings.timesteps):    
            if self.sim_a_0 <= 0:                         # All reactants got exhausted
                settings.endtime = 10**50
                break
            
            self.RunExactTimestep()                       # Run direct SSA 
            self.HandleEvents()
            
            # Update Propensities selectively            
            if self.sim_t < settings.endtime:  
                if not self._IsPerformEvent:
                    self.species_to_update = self.parse.reaction_affects[self.reaction_index] # Determine vars to update                
                else:
                    self.species_to_update = [s for s in range(self.n_species)]         
                
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
            sys.stdout.write('\rsimulation done!               \n')
            

    def RunExactTimestep(self):
        """ Calculates a time step of the Direct Method """ 
        if self.sim_t == 0:
            randoms = np.random.random(1000) 
            self.randoms_log = np.log(randoms)*-1
            self.randoms = np.random.random(1000)
            self.count = 0           
        elif self.count == 1000:
            randoms = np.random.random(1000) 
            self.randoms_log = np.log(randoms)*-1
            self.randoms = np.random.random(1000)    
            self.count = 0      
        
        self.sim_tau = self.randoms_log[self.count]/float(self.sim_a_0)       # reaction time generation
        self.sim_r2 = self.randoms[self.count]                                # Draw random number 2 [0-1]
        self.count +=1
        
        if (self.sim_t + self.sim_tau) < self.settings.endtime:
            self.sim_t += self.sim_tau                                        # Time update
            self.reaction_index = 0
            sum_of_as = self.sim_a_mu[self.reaction_index]
            criteria = self.sim_r2*self.sim_a_0
            while sum_of_as < criteria:                                       # Use r2 to determine which reaction will occur
                self.reaction_index += 1	                                  # Index
                sum_of_as += self.sim_a_mu[self.reaction_index]                 

            try:
                self.X_matrix += self.N_matrix_transpose[self.reaction_index]
                self.timestep += 1
            except MemoryError as ex:
                print(ex)
                sys.exit()               
        else:            
            self.sim_t = self.settings.endtime            
            self.reaction_index = np.nan
