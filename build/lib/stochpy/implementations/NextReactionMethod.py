#! /usr/bin/env python
"""
Next Reaction Method
====================
This module performs the Next Reaction Method from Gibson and Bruck [1]. Therefore, it is also called the Gibson and Bruck algorithm.

[1] M.A. Gibson and J. "Bruck Efficient Exact Stochastic Simulation of Chemical Systems with Many Species and Many
Channels", J. Phys. Chem., 2000, 104, 1876-1889

Written by T.R. Maarleveld, Amsterdam, The Netherlands, Speed improvements by M. Moinat.
E-mail: tmd200@users.sourceforge.net
Last Change: June 23, 2015
"""

############################# IMPORTS #################################### 

from __future__ import division, print_function, absolute_import

import time,os,sys,copy

from ..tools.Priority_queue import PQDict
from .StochPyTools import  __species__,StochPySSA_Shared,np

########################### END IMPORTS ##################################

class NextReactionMethod(StochPySSA_Shared):
    """ 
    Next Reaction Method from Gibson and Bruck 2000 [1]. 

    [1] M.A. Gibson and J. "Bruck Efficient Exact Stochastic Simulation of Chemical Systems with Many Species and Many Channels", J. Phys. Chem., 2000,104,1876-1889 

    Input:  
     - *model_file* filename.psc
     - *model_dir* /home/user/Stochpy/pscmodels/filename.psc    
    """
    def __init__(self,model_file,model_dir,IsQuiet=False):
        self.Parse(model_file,model_dir,IsQuiet=IsQuiet,IsNRM=True,IsDelayed = False, IsSMM = False)  
              

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
        self.X_matrix = copy.deepcopy(settings.X_matrix)  
        self.fixed_species_amount = copy.deepcopy(self.parse.fixed_species_amount)
       
        try:
            self.volume_code = settings.volume_code
        except AttributeError: 
            self.volume_code = "self._current_volume = 1"    # No volume_code present in settings (normal mode)        
        
        self.BuildPropensityCodes()          
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
            self.InitialMonteCarloStep()    # Initialize Tau Heap 
            if not settings.IsOnlyLastTimepoint:
                self.Initial_Conditions()
        else:
             #print(self._IsPerformEvent)             
             self.UpdateHeap()              # Update the heap because of the changed propensities after e.g. division            

        nstep_counter = 1
        t1 = time.time()
        while self.sim_t < self.settings.endtime and self.timestep < self.settings.timesteps:            
            if self.sim_a_0 <= 0:           # All reactants got exhausted
                 settings.endtime = 10**50
                 break            
                  
            self.RunExactTimestep()         # Run time step   
            self.HandleEvents()               
            
            if self.sim_t < settings.endtime: 
                self.Propensities()         # Update Propensities  
                self.UpdateHeap()                                           
                
                if not settings.IsOnlyLastTimepoint:      
                    self.GenerateOutput()   # Store Output    
                    
            self._IsPerformEvent = False    # set to false (or just to make sure).            
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
            
   
    def InitialMonteCarloStep(self): 
        """ Monte Carlo step to determine all taus and to create a pqdict indexed priority queue """  
        randoms = np.random.random(self.n_reactions)                     # Pre-generate for each reaction random numbers
        self.randoms_log_init = -1 * np.log(randoms)      
        self.sim_taus = self.randoms_log_init/self.sim_a_mu + self.sim_t # Make all taus absolute with sim_t (= given starttime)
        
        pairs = [(j, self.sim_taus[j]) for j in range(self.n_reactions)] #(key, priority) pairs   
        self.heap = PQDict(pairs)        
        
        ##Build initial randoms for tau's to come, always at start of execution NRM
        self.randoms_log = -1 * np.log(  np.random.random(1000)  )       # Pre-generate random numbers
        self.count = 0
        
        
    def UpdateHeap(self): #03-11-2013
        """ Renews the changed propensities in the priority queue heap. """
        if not self._IsPerformEvent:
            self.prop_to_update = self.parse.dep_graph[self.reaction_index] # Propensities to update  
        else:
            self.prop_to_update = [r for r in range(self.n_reactions)]        
              
        if self.count >= (1000-len(self.prop_to_update)):
            randoms = np.random.random(1000)                             # Pre-generate random numbers
            self.randoms_log = -1 * np.log(randoms)
            self.count = 0
            
        for j in self.prop_to_update: #Updated propensities should also be updated in the heap
            if j == self.reaction_index: #5c
                tau_new = self.randoms_log[self.count]/self.sim_a_mu[j] + self.sim_t
                self.count += 1
            else: #5b, changed propensity
                if self.sim_a_mu_prev[j] == 0: 
                    ##Note: due to zero propensity (and an inf tau), it is difficult to reuse random. Assume getting new random is faster.
                    tau_new = self.randoms_log[self.count]/self.sim_a_mu[j] + self.sim_t
                    self.count += 1
                else: #Normal update
                    tau_alpha = self.sim_taus[j]        #Faster than getting from heap directly (self.heap[j])
                    tau_new = self.sim_a_mu_prev[j]/self.sim_a_mu[j]*(tau_alpha- self.sim_t) + self.sim_t              
            #Note, no exception for self.sim_a_mu == 0. Then tau_new automatically becomes infinite (faster than if statement to except this)          
            self.sim_taus[j] = tau_new #Keep track of the reaction taus, parallel to the heap.
            self.heap.updateitem(j, tau_new) #Updates the tau with index j and resorts the heap.
            

    def Propensities(self):
        """ Determines the propensities to fire for each reaction at the current time point. At t=0, all the rate equations are compiled. """
        #print(self._IsInitial)
        if self._IsInitial:
            self.sim_a_mu = np.zeros(self.n_reactions)            # Initialize a(mu)
            [setattr(__species__,self.parse.species[s],self.X_matrix[s]) for s in range(self.n_species)] # Set species quantities
            [setattr(__species__,self.fixed_species[s],self.fixed_species_amount[s]) for s in range(len(self.fixed_species))]
            self.reaction_fired = -1                              # Update all propensities     
            self._IsInitial=False           
        else:          
            self.sim_a_mu_prev = copy.copy(self.sim_a_mu)         # Backup old propensity 
            if self._IsPerformEvent:
                [setattr(__species__,self.parse.species[s],self.X_matrix[s]) for s in range(self.n_species)] #Update all species, to be sure.
                self.reaction_fired = -1                          # Update all propensities                
            else:
                self.species_to_update = self.parse.reaction_affects[self.reaction_index] # Determine vars to update
                [setattr(__species__,self.parse.species[s],self.X_matrix[s]) for s in self.species_to_update]
                self.reaction_fired = self.reaction_index         # Update the propensities that depend on this reaction
        
        propensity_eval_code = self.propensity_codes[ self.reaction_fired ] #21-11-2013, select code of subset to be updated. [-1] updates all
        self.rateFunc(propensity_eval_code, self.sim_a_mu)        # Calc. Propensities and put result in sim_a_mu
        #print(self.X_matrix,self.sim_a_mu)
        assert self.sim_a_mu.min() >= 0, "Error: Negative propensities are found" 
        self.sim_a_mu = abs(self.sim_a_mu)                        # -0 to 0
        self.sim_a_0 = self.sim_a_mu.sum()
        

    def RunExactTimestep(self):
        """ Perform a direct SSA time step and pre-generate M random numbers """        
        self.minimum = self.heap.peek()                           # peek() returns item at top of heap (has lowest tau)            
        self.reaction_index = self.minimum[0]                     # Pick reaction to executeO(1)
        self.sim_tau = self.minimum[1]                            # Pick tau O(1)
        
        if self.sim_tau < self.settings.endtime:
            self.sim_t = self.sim_tau                             # New time
            try:
                self.X_matrix += self.N_matrix_transpose[self.reaction_index]
                self.timestep += 1
            except MemoryError as ex:
                print(ex)
                sys.exit()            
        else:           
            self.sim_t = self.settings.endtime
            self.reaction_index = np.nan            
