#! /usr/bin/env python
"""
Delayed Direct Method
=====================

This module performs the Delayed Direct Stochastic Simulation Algorithm adapted from Cai (2007) [1].

This algorithm is used to generate exact realizations of the Markov jump process with delays. Of course, the algorithm is stochastic, so these realizations are different for each run.

Only molecule populations are specified. Positions and velocities, such as in Molecular Dynamics (MD) are ignored. This makes the algorithm much faster, because non-reactive molecular collisions can be ignored.different
Still, this exact SSA is quite slow, because it insists on simulating every individual reaction event, which takes a lot of time if the reactant population is large. Furthermore, even larger problems arise if the model contains distinct processes operating on different time scales [2].

[1] Xiaodong Cai (2007), "Exact stochastic simulation of coupled chemical reactions with delays", J.Phys. Chem. 126:124108

Written by M Moinat, T.R. Maarleveld
Adapted from 'DirectMethod' by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: June 01, 2015
"""

from __future__ import division, print_function, absolute_import

__doc__ = """
          Delayed Direct Method
          =====================
          This module performs the Delayed Direct Stochastic Simulation Algorithm adapted from Cai (2007) [1]. This algorithm is used to generate exact realizations of the Markov jump process with delays. Of course, the algorithm is stochastic, so these realizations are different for each run.
          Only molecule populations are specified. Positions and velocities, such as in Molecular Dynamics (MD) are ignored. This makes the algorithm much faster, because non-reactive molecular collisions can be ignored.
          Still, this exact SSA is quite slow, because it insists on simulating every individual reaction event, which takes a lot of time if the reactant population is large. Furthermore, even larger problems arise if the model contains distinct processes operating on different time-scales [2].
          
          [1] Xiaodong Cai (2007), "Exact stochastic simulation of coupled chemical reactions with delays", J.Phys. Chem. 126:124108
          """
############################# IMPORTS ####################################

import sys,copy,time,os,operator
from .StochPyTools import __species__,StochPySSA_Shared,np
    
########################### END IMPORTS ##################################

class DelayedDirectMethod(StochPySSA_Shared):
    """ 
    Delayed Direct Stochastic Simulation Algorithm from Cai (2007) [1].

    This algorithm is used to generate exact realizations of the Markov jump process with delays. Of course, the algorithm is stochastic, so these realizations are different for each run.

    [1] Xiaodong Cai (2007), "Exact stochastic simulation of coupled chemical reactions with delays", J.Phys. Chem. 126:124108

    Input:  
     - *File* filename.psc
     - *dir* /home/user/Stochpy/pscmodels/filename.psc    
    """
    def __init__(self,File,dir,IsQuiet=False):        
        self.Parse(File,dir,IsQuiet=IsQuiet,IsNRM = False, IsDelayed=True, IsSMM = True)     
        
 
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
        except AttributeError:            # No volume_code present in settings, Volume always 1
            self.volume_code = "self._current_volume = 1"
        
        self._IsTimeEvent = False        
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
            # Delayed Direct method specific         
            self.Tstruct = [(0, np.nan), (np.inf, np.nan)]  # Will hold delayed reactions.
            self._Tstruct_prev =  [(0, np.nan), (np.inf, np.nan)]              

        nstep_counter = 1
        t1 = time.time()
        while (self.sim_t < settings.endtime) and (self.timestep < settings.timesteps):          
            # If all propensities zero and no delayed reactions pending.            
            if self.sim_a_0 == 0 and len(self.Tstruct) <= 2:   
                settings.endtime = 10**50             
                break

            self.RunExactTimestep()                                          # Run direct SSA 
            self.HandleEvents()
            
            # Update Propensities selectively
            if self.sim_t < settings.endtime:
                if self._IsPerformEvent:
                    self.species_to_update  = [s for s in range(self.n_species)]                  
                    if self._IsTimeEvent:                                    # Reset the Tstruct to before event trigger, added delay will be removed. Only happens for time events, for copy number events the delay should be kept in the Tstruct.
                        self.Tstruct =  copy.copy(self._Tstruct_prev)
                        # Update relative delay time to new sim_tau (time up to event)
                        true_sim_tau = self.sim_t - self._current_sim_t     #10-04-2014
                        self.Tstruct = [self.Tstruct[0]] + [(T - true_sim_tau, j) for (T,j) in self.Tstruct[self.i+1:]] # 11-04-2014 self.i to remove finished delayed reactions.
                        self._IsTimeEvent = False        
                self.Propensities()
            
                if not settings.IsOnlyLastTimepoint:
                    self.GenerateOutput()
                
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
        #np.random.seed(5)
        if self.sim_t == 0:
            self.randoms1 = np.random.random(1000)
            self.randoms2 = np.random.random(1000)
            self.count = 0       
        elif self.count == 1000:
            self.randoms1 = np.random.random(1000)
            self.randoms2 = np.random.random(1000)
            self.count = 0       
            
        self.GenerateTau()                                                      # reaction time generation from self.randoms1
        if not self._IsPerformEvent:                                            # If event triggered by delay, this replaces the initiation of reaction. 
            self.sim_t += self.sim_tau           
            if self.sim_a_0 == 0:                                               # No reaction initiates (but delayed reaction can still complete in next iteration).
                self.reaction_index = np.nan                
            elif self.sim_t < self.settings.endtime:
                self.InitiateReaction()                                         # Pick reaction and perform it
                #TODO: reaction (sim_t) after event, reaction will never occur but delayed event is added. Problem: Delayed reaction started, at time of event..?
                self.count += 1
                self.timestep += 1            
            else: 
                self.sim_t = self.settings.endtime
                self.reaction_index = np.nan
                
            
    def GenerateTau(self): 
        """ Generates the waiting time, according to the improved equation 11 of Cai """
        self.sim_r1 = self.randoms1[self.count]     
        self._current_sim_t = self.sim_t                                        # Event handling: Remember sim_t at current timepoint.        
        self.i = 0
        if len(self.Tstruct) == 2:                                              # No ongoing delayed reaction
            self.sim_tau = -1 * np.log(self.sim_r1)/float(self.sim_a_0)           
        else:            
            at_previous = 0
            self.at = self.sim_a_0 * self.Tstruct[self.i+1][0]                  # Up to first delayed reaction
            F = 1 - np.exp(-self.at)
            while F < self.sim_r1 and not self._IsPerformEvent:                 # End if initiation of reaction reached, or if event started.
                self.i += 1
                #Stop while loop if over end_time
                if self.sim_t + self.Tstruct[self.i][0] >= self.settings.endtime:
                    if self.Tstruct[self.i][0] == np.inf: #09-01-2014
                        self.sim_tau = self.Tstruct[self.i - 1][0]              # Time to last delayed reaction, prevents entime of 1e10 if mode is steps.
                    else:
                        self.sim_tau = self.settings.endtime - self.sim_t       # Time to endtime
                    self.i -= 1                                                 # The ith delayed reaction has never finished (is after end_time)
                    break 
                    
                at_previous = self.at
                #Update species and propensities (also a_0) due to completion of delayed reaction at Ti
                self.reaction_index = self.Tstruct[self.i][1]
                self.sim_t += self.Tstruct[self.i][0]                           # 11-2-2014 Temporary sim_t update                   

                ## Handle time events initiated by delayed reaction. This event is executed instead of delay, but delayed reaction will be stored.
                self.HandleEvents()                                             # 28-03-2014
                if self._IsPerformEvent:                    
                    self.sim_tau = self.sim_t - self._current_sim_t             # Time difference last event, needed for update Tstruct
                    self.i -= 1                   
                    break
                    
                self.CompleteReaction()                                         # TODO: copy number initiation due to completion delayed. Again check for event after this?
                self.Propensities()  
                
                ## Handle copy number events initiated by completion of delayed reaction. 
                self.HandleEvents()                                             # 28-03-2014
                if self._IsPerformEvent:                    
                    self.sim_tau = self.sim_t - self._current_sim_t             # Time difference last event, needed for update Tstruct
                    break
                    
                if not self.settings.IsOnlyLastTimepoint:                       # Only save output if not specified that only last output should be saved.
                    self.GenerateOutput( completion_delayed = True )            # timepoint is the absolute time of completion                    
                self.sim_t -= self.Tstruct[self.i][0]
                self.timestep += 1
                    
                #Test whether reaction time (self.sim_r1) falls after next delayed reaction (If not, the while loop will fail)
                if self.Tstruct[self.i+1][0] != np.inf:
                    self.at += self.sim_a_0 * (self.Tstruct[self.i+1][0] - self.Tstruct[self.i][0])
                    F = 1 - np.exp(-self.at)
                else:
                    F = np.inf
                
                if self.reactions_NonConsuming:                                 # Saves time
                    self.CheckReactantExhaustion()                              # Deletes pending reactions from Tstruct if reaction exhausted.
                
            else:                                                               # Only execute when F<self.sim_r1, not if break encountered
                if self.sim_a_0 == 0:                                           # This will be true if no reactants AND no delayed reactions pending.
                    self.sim_tau = self.Tstruct[self.i][0]                      # The simulation has to end after the last delayed completion
                else:
                    self.sim_tau = self.Tstruct[self.i][0] - ( np.log(1-self.sim_r1) + at_previous )/float(self.sim_a_0)
            
            # Store Tstruct before modified. Needed to revert new delayed reaction and time update of delayed reactions if event is triggered.
            self._Tstruct_prev  =  copy.copy(self.Tstruct)                      #10-04-2014
            #Remove completed (all upto place self.i+1) and update times to completion (bring forward)
            self.Tstruct = [self.Tstruct[0]] + [(T - self.sim_tau, j) for (T,j) in self.Tstruct[self.i+1::]]  
               
                
    def InitiateReaction(self):
        """ Chooses reaction that will occur at sim_t and updates species accordingly """
        self.sim_r2 = self.randoms2[self.count]                                 # Draw random number 2 [0-1]
        self.reaction_index = 0
        sum_of_as = self.sim_a_mu[self.reaction_index]
        criterium = self.sim_r2*self.sim_a_0
        while sum_of_as < criterium:                                            # Use r2 to determine which reaction will occur
            self.reaction_index += 1                                            # Index
            sum_of_as += self.sim_a_mu[self.reaction_index]    
        
        # Update species according to consuming/nonconsuming/normal
        if self.reaction_index in self.reactions_Consuming:                     # Only update reactants
            self.X_matrix += self.N_matrix_transpose_reactants[self.reaction_index]      
            self.species_to_update = self.parse.reactant_indices[ self.reaction_index ]
            self.add_delay()
            
        elif self.reaction_index in self.reactions_NonConsuming:                      
            self.add_delay()                                                    # After the delay there will be completion of consumption and production.
            self.species_to_update = []                                         # No species have to be updated, as there is no change in species upon initiation of nonconsuming.
        
        else:                                                                   # A normal, nondelayed reaction.
            self.X_matrix += self.N_matrix_transpose[self.reaction_index] 
            self.species_to_update = self.parse.reaction_affects[ self.reaction_index ]
            

    def CompleteReaction(self):
        """ Completes a delayed reaction and updates propensities accordingly."""        
        if self.reaction_index in self.reactions_Consuming:                    
            self.X_matrix += self.N_matrix_transpose_products[self.reaction_index]
            self.species_to_update = self.parse.product_indices[ self.reaction_index ]# Only complete production (reactants already consumed)            
        elif self.reaction_index in self.reactions_NonConsuming: 
            self.X_matrix += self.N_matrix_transpose[self.reaction_index]  
            self.species_to_update = self.parse.reaction_affects[ self.reaction_index ]        # Both reactants and products complete at the same time            
        else:
            print("Error: This should not happen; a non-delayed reaction is completing with a delay...")
                
    def add_delay(self):
        """ Add delay """
        #Delay from distribution specified in DoDelayedStochSim (StochSim.py)
        delay = self.distr_functions[ self.reaction_index ]( *self.distr_parameters[ self.reaction_index ] ) #Each delayed reaction has its own distribution        
        if type(delay) == np.ndarray:
            print("Warning: Delay parameters do not match the distribution ({0} with {1})".format(self.distr_functions[ self.reaction_index ], self.distr_parameters[ self.reaction_index ]))
            delay = delay[0]
        
        #Check for negative delays
        i=0
        while delay < 0:
            print("Warning: The chosen distribution ({0}) produced a negative delay. Drawing new delay.".format(self.distr_functions[ self.reaction_index ]) )
            if i == 0:
                print("Info: Delay distributions will be distorted and simulation speed will be drastically lower.")
            delay = self.distr_functions[ self.reaction_index ]( *self.distr_parameters[ self.reaction_index ] )
            i+=1
            if i >= 10:
                print("Error: negatives delays keep occurring. Please choose different delay distributions.")
                sys.exit(1)
        
        # Insert delay and resort
        self.Tstruct.insert(-1, (delay, self.reaction_index) )          
        self.Tstruct.sort()       
    
    def CheckReactantExhaustion(self):
        """ Checks wheter a reactant is exhausted. Important if still reactions pending in Tstruct, but no reactant present. """
        for r_nonconsuming in self.reactions_NonConsuming:
            for r in self.parse.reactant_indices[r_nonconsuming]:
                if self.X_matrix[r] == 0:                                       # Reactant Exhausted
                    # Remove all pending r_nonconsuming reactions from Tstruct. Order is not affected.
                    temp_Tstruct1 = copy.copy(self.Tstruct[:self.i+1])          # Don't alter processed part of Tstruct
                    temp_Tstruct2 = copy.copy(self.Tstruct[self.i+1:])          # Only remove from non-processd Tstruct
                    temp_Tstruct2 = [tau_pair for tau_pair in temp_Tstruct2 if tau_pair[1] != r_nonconsuming] #list(filter(lambda tau_pair: tau_pair[1] != r_nonconsuming, temp_Tstruct2))  Python 3.x                   
                    self.Tstruct = temp_Tstruct1 + temp_Tstruct2  
