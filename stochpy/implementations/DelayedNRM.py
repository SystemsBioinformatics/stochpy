#! /usr/bin/env python
"""
Modified Next Reaction Method with delays
=========================================

This module performs the Next Reaction Method with delay from Anderson [1], Algorithm 7.
Note that the definition of Cai is used for types of delayed reactions: consuming and non-consuming.

[1] David F. Anderson "A modified next reaction mehtod for simulating chemical systems with 
time dependent propensities and delays", J. Phys. Chem., 2007, 127, 214107 

Author: M. Moinat, T.R. Maarleveld
Adapted from 'NextReactionMethod' by T.R. Maarleveld, Amsterdam, The Netherlands
Last Change: June 08, 2015
"""
from __future__ import division, print_function, absolute_import
import sys,copy,time,os

from .StochPyTools import __species__,StochPySSA_Shared,np

#############################################################

class DelayedNRM(StochPySSA_Shared):
    """ 
    Modified Next Reaction Method from Anderson 2007 [1]. 

    [1] David F. Anderson "A modified next reaction mehtod for simulating chemical systems with 
      time dependent propensities and delays", J. Phys. Chem., 2007, 127, 214107

    Input:  
     - *model_file* filename.psc
     - *model_dir* /home/user/Stochpy/pscmodels/filename.psc    
    """
    def __init__(self,model_file,model_dir,IsQuiet=False):
        self.Parse(model_file,model_dir,IsQuiet=IsQuiet,IsNRM = True, IsDelayed=True,IsSMM = False)   
           
    def Execute(self,settings,IsStatusBar=False):
        """
        Generates steps of the Markov jump process.

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
        except AttributeError:                             # No volume_code present in settings
            self.volume_code = "self._current_volume = 1"
        
        #self.species_to_update = [s for s in range(self.n_species)] # ensure that the first run updates all species   
        self.BuildInits()                                  # Step 1        
        self.Propensities()                                # Step 2      

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
            # Delayed Direct method specific         
            self.Tstruct = [(0, np.nan), (np.inf, np.nan)]  # Will hold delayed reactions      
            
        #Step 3 Initialize P        
        m_randoms = np.random.random(size = self.n_reactions)
        self.P = -1*np.log(m_randoms)                       # Exponential (unit rate) distribution
        self.T = np.zeros(self.n_reactions)                 # Step 1
        self.BuildTaus()                                    # Step 4, initialization        
        
        nstep_counter = 1
        t1 = time.time()
        while self.sim_t < settings.endtime and self.timestep < settings.timesteps:
            if self.sim_a_mu.sum() == 0 and len(self.Tstruct) <= 2:       # All reactants got exhausted and no delayed reactions pending.
                 settings.endtime = 10**50
                 break      

            self.RunExactTimestep()                         # Step 5-10             
            self.HandleEvents()

            #Update system, only if max time or steps is not passed
            if (self.sim_t < settings.endtime) and (self.timestep < settings.timesteps):                
                # Bug Fix Maxim september 05, 2014
                if self._IsPerformEvent:             
                    m_randoms = np.random.random(size = self.n_reactions)
                    self.P = -1*np.log(m_randoms)           # Exponential (unit rate) distribution
                    self.T = np.zeros(self.n_reactions)     # Step 1
                    self.BuildTaus()                        # Step 4, initialization        
                    self.species_to_update  = [s for s in range(self.n_species)]   
                else:
                    self.UpdateTP()                         # Steps 11,12                              
                
                self.Propensities()                         # Step 13
                
                #Update taus (and order)
                self.BuildTaus()                            # Step 4, iterative                        
                  
                if not settings.IsOnlyLastTimepoint: #16-01-2014
                   self.GenerateOutput()                
                    
            if (self.sim_t < settings.endtime):         
                self.timestep += 1
            t2 = time.time()
            if IsStatusBar and t2-t1> 1:
                t1 = time.time()
                sys.stdout.write('\rsimulating {0:s}\r'.format('.'*nstep_counter) ) 
                sys.stdout.flush() 
                nstep_counter+=1
                if nstep_counter > 10:
                    nstep_counter = 1
                    
        if settings.IsOnlyLastTimepoint or settings.endtime != 10**50: 
            self.GenerateOutput()
                    
        if IsStatusBar and t1 and not settings.quiet:  
            sys.stdout.write('\rsimulation done!               \n')                          
          
    def BuildInits(self):
        """ Build initials that are necessary to generate a trajectory """    
        self.rFire= np.ones(self.n_reactions)
        self.sim_a_mu_zero = []
        for i in range(self.n_reactions):
            self.sim_a_mu_zero.append({'t1':0,'a_alpha_old':0,'tau_alpha':0})    
        self.sim_taus = None                            # necessary for iteration

    def BuildTaus(self):
        """ (Re)Calculates and orders taus """
        self.sim_taus = (self.P - self.T)/self.sim_a_mu # Calculates all initial taus from P and T arrays.
        self.tauPairs = [(tau, reaction_index) for reaction_index, tau in enumerate(self.sim_taus)]      
        self.tauPairs.sort()                            # Sorts on reaction times (taus), ascending.
    
    def UpdateTP(self):
        """ Update T and P """
        #Update all T's (= a kind of integral)
        self.T += self.sim_a_mu * self.tau              # Step 8
        
        #Only update P of reaction that initiated (so is not a completed delayed reaction). Only random generated each step is here
        if not self._IsCompletionDelayed:
            random = np.random.random()     
            self.P[self.reaction_index] += -1*np.log(random)
            
    def RunExactTimestep(self):
        """ Perform a direct SSA time step """     
        self.GetTauReaction()                           # Step 5        
        self.sim_t += self.tau                          # Step 6            
        if self.sim_t < self.settings.endtime:
            #Perform reaction
            if self._IsCompletionDelayed:               # Step 7, complete the delayed reaction.
                self.CompleteReaction()         
                if self.reactions_NonConsuming:         # Saves time if no nonconsuming
                    self.CheckReactantExhaustion()      # 9-01-2014 Prevents negative species due to delayed nonconsuming reactions.
            else:
                self.InitiateReaction()                 # Steps 8,9,10

            self.prop_to_update = self.parse.dep_graph[self.reaction_index] # Propensities to update, does not see difference N_matrix_reactants or N_matrix_products into account      
            if self.reaction_index not in self.prop_to_update: 
                self.prop_to_update.append(self.reaction_index)     # Always add fired reaction to update  
            self.rFire= np.zeros(self.n_reactions)    
            for pos in self.prop_to_update: 
                self.rFire[pos] = 1                                 # Makes sure that in Propensities() only changed propensities get updated
        else: 
            self.sim_t = self.settings.endtime
            self.reaction_index = np.nan
            
    def GetTauReaction(self):
        """ Gets minimum tau, corresponding reaction index and whether it initiates or completes """
        #Tau is relative, not absolute as in normal NRM (!!!)
        self.minimum_reaction = self.tauPairs[0]                   # Pick first initiation reaction to occur                        
        self.tau_reaction = self.minimum_reaction[0]               # Pick tau of the reaction 
        
        self.minimum_delay = self.Tstruct[1]                       # From sorted Tstruct (on index 0 is the initial (0, np.nan))
        self.tau_delay = self.minimum_delay[0] - self.sim_t        # Convert from absolute to relative time            
        
        if self.tau_reaction <= self.tau_delay:                    # Initiate a reaction
            self.tau = self.tau_reaction
            self.reaction_index = self.minimum_reaction[1]
            self._IsCompletionDelayed = False
        else:                                                      # A delayed reaction completes (self.tau_delay < self.tau_reaction)
            self.tau = self.tau_delay
            self.reaction_index = self.minimum_delay[1]
            self._IsCompletionDelayed = True      
    
    def InitiateReaction(self):
        """ Initiates a reaction based on being delayed (consuming or nonconsuming) or not """
        if self.reaction_index in self.reactions_Consuming:                     # Step 10, Consuming delayed
            self.X_matrix += self.N_matrix_transpose_reactants[self.reaction_index]             # Updates only reactants
            self.species_to_update = self.parse.reactant_indices[ self.reaction_index ]
            self.add_delay()
        elif self.reaction_index in self.reactions_NonConsuming:                # Step 9, Non-consuming delayed         
            self.add_delay()
            self.species_to_update = []   
        else:                                                                   # Step 8, A normal, non-delayed reaction
            self.X_matrix += self.N_matrix_transpose[self.reaction_index]
            self.species_to_update = self.parse.reaction_affects[ self.reaction_index ]
    
    def CompleteReaction(self):
        """ Completes a delayed reaction, depending on consuming or not or not"""      
        if self.reaction_index in self.reactions_Consuming:                     # Only updates products
            self.X_matrix += self.N_matrix_transpose_products[self.reaction_index]
            self.species_to_update = self.parse.product_indices[ self.reaction_index ]            
        elif self.reaction_index in self.reactions_NonConsuming:                # Updates reactants and products
            self.X_matrix += self.N_matrix_transpose[self.reaction_index]
            self.species_to_update = self.parse.reaction_affects[ self.reaction_index ]            
        else:
            print("This should not happen. A non-delayed reaction is completing with a delay...")
            
        self.Tstruct.pop(1)                                                     # Remove the completion time from the list, assumes Tstruct is sorted and Tstruct[0] = (0,np.nan)
        
    def add_delay(self):
        """ Adds delay time from distribution, specified in DoDelayedStochSim (StochSim.py) """
        delay = self.distr_functions[ self.reaction_index ]( *self.distr_parameters[ self.reaction_index ] )    # Each delayed reaction has its own distribution

        if type(delay) == np.ndarray:
            print("Warning! Delay parameters don't match the distribution ({0} with {1})".format(self.distr_functions[ self.reaction_index ], self.distr_parameters[ self.reaction_index ]))
            delay = delay[0]
        
        #Check for negative delays
        i=0
        while delay < 0:
            print("Warning: The chosen distribution ({0}) produced a negative delay. Drawing new delay.".format(self.distr_functions[ self.reaction_index ]))
            if i==0:
                print("Delay distributions will be distorted and simulation speed will be drastically lower.")
            delay = self.distr_functions[ self.reaction_index ]( *self.distr_parameters[ self.reaction_index ] )
            i += 1
            if i >= 10:
                print("Error: negatives delays keep occuring. Please choose different delay distributions.")
                sys.exit(1)
        
        #Add ABSOLUTE completion time and resort Tstruct
        self.Tstruct.insert(-1, (delay + self.sim_t, self.reaction_index)) 
        self.Tstruct.sort()
    
    def CheckReactantExhaustion(self):
        """ Checks wheter a reactant is exhausted after completion and deletes the pending delayed reactions of this reactant."""
        for r_nonconsuming in self.reactions_NonConsuming:
            for r in self.parse.reactant_indices[r_nonconsuming]:
                if self.X_matrix[r] == 0:                     # Reactant Exhausted
                    # Remove all pending r_nonconsuming reactions from Tstruct, if one of reactants exhausted. Order is not affected.
                    self.Tstruct = [tau_pair for tau_pair in self.Tstruct if tau_pair[1] != r_nonconsuming]   
                                   #list(filter(lambda tau_pair: tau_pair[1] != r_nonconsuming, self.Tstruct)) Python 3.x   
