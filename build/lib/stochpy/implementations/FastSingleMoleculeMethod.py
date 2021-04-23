"""
Fast Single Molecule Method (fSMM)
==================================

This module performs the fast Single Molecule Method, based on the Next Reaction Method by Gibson and Bruck [1]. 
Only user specified reaction will be treated as single molecule reactions, others will be treated as 'normal' exponential reactions. 
This method is called the fSMM since only zero and first-order single molecule reactions can be used. 'normal' exponential reactions can be from any order.

[1] M.A. Gibson and J. "Bruck Efficient Exact Stochastic Simulation of Chemical Systems with Many Species and Many
Channels", J. Phys. Chem., 2000, 104, 1876-1889

Note: do not use net stoichiometries

Written by M Moinat, T.R. Maarleveld
Adapted from NextReactionMethod.py by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: June 08, 2015
"""

############################## IMPORTS ###################################

from __future__ import print_function

import sys,copy,time,os
from ..tools.Priority_queue import PQDict
from .StochPyTools import __species__,StochPySSA_Shared,np

########################### END IMPORTS ##################################

class FastSingleMoleculeMethod(StochPySSA_Shared):
    """
    This module performs a Single Molecule Method, based on the Next Reaction Method by Gibson and Bruck [1]. 

    [1] M.A. Gibson and J. "Bruck Efficient Exact Stochastic Simulation of Chemical Systems with Many Species and Many Channels", J. Phys. Chem., 2000,104,1876-1889 

    Input:  
     - *model_file* filename.psc
     - *model_dir* /home/user/Stochpy/pscmodels/filename.psc    
    """
    def __init__(self,model_file,model_dir,IsQuiet=False):
        self.Parse(model_file,model_dir,IsQuiet=IsQuiet,IsTauleaping = False, IsNRM = True, IsDelayed=False, IsSMM=True)      
    
    def Execute(self,settings,IsStatusBar=False):
        """
        Generates T trajectories of the Markov jump process.

        Input:
         - *settings* (class object)   
        """
        if settings.IsSeed:
            np.random.seed(5)
        
        self._IsInitial = True
        self.sim_t = copy.copy(settings.starttime)
        self.X_matrix = copy.deepcopy(settings.X_matrix)
        self.fixed_species_amount = copy.deepcopy(self.parse.fixed_species_amount)
            
        self.settings = copy.copy(settings)                
        
        self._IsExhausted = False # 10-01-2014 moved from runexacttimestep
        self.single_molecule_reactions = list(self.distr_functions)        # Reactions with assigned distribution function (these are molecular)
        
        try:
            self.volume_code = settings.volume_code
        except AttributeError: 
            self.volume_code = "self._current_volume = 1"    # No volume_code present in settings (normal mode)        

        self.BuildPropensityCodes()                       
        self.BuildDependencyGraphs()        
        self.Propensities()
        self.InitializeTauHeap()           

        if not self.sim_t:
            self.timestep = 1   
            self.sim_output = []
            self.propensities_output = []              
            self.V_output = []        
            self._IsTrackPropensities = False # TODO: check      
            self.SpeciesSelection() 
            #self.RateSelection()
            self.SetEvents()                 
            if not settings.IsOnlyLastTimepoint:
                self.Initial_Conditions()      
            
        nstep_counter = 1
        t1 = time.time()
        while self.sim_t < settings.endtime and self.timestep < settings.timesteps:              
            if self._IsExhausted:                          # All reactants got exhausted
                settings.endtime = 10**50
                break

            self.RunExactTimestep(settings)                # Run time step     
            self.HandleEvents()     
            
            if not self._IsPerformEvent:                   # 28-03-2014
                if (self.sim_t < settings.endtime) and (self.timestep < settings.timesteps):     
                    self.UpdateSystem()
            else:
                self.Propensities()
                self.InitializeTauHeap()
            
            if not settings.IsOnlyLastTimepoint:
                self.GenerateOutput()  
            
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
        
    def BuildDependencyGraphs(self): 
        """
        Three different dependency graphs exist: G, H, and I of which the latter two are dedicated to the single molecule reactions (SMRs).
        - G is the same dependency graph as used in the Next Reaction Method
        - H (reactant, reactant): contains for every reaction k the SMRs that are affected by the loss of reactants due to firing of reaction k (reaction k itself is excluded)
        - I (product, reactant): contains for every reaction k the SMRs that are affected by the formation of products (i.e. reactants for other SMRs) due to firing of reaction k.        
        
        Note: For the formation of H we have to use "reactants" rather than "depends_on". The latter contains also species that catalyze the reaction, so they are not lost after firing of the reaction. This means that we do not have to delete putative waiting times for reactions where this species acts as reactant.
        """          
        # H: Intersection reactants[i] and reactants[j]. If reaction i fires, which SMR j needs to have waitingtime deleted.
        self.dep_graph_H = []
        for i,reactants_i in enumerate(self.parse.reactant_indices):
            to_change = []
            for j,reactants_j in enumerate(self.parse.reactant_indices):
                if j not in self.single_molecule_reactions or i == j:   # Only include SMRs and not own index (processed separately)
                    continue
                for reactanti in reactants_i:                    
                    if reactanti in reactants_j:
                        to_change.append(j)
            self.dep_graph_H.append(to_change)
        
        # I: Intersection products[i] and reactants[j]. If reaction i fires, which single molecule reaction j needs to have waitingtime added.
        self.dep_graph_I = []
        for i,reactants_i in enumerate(self.parse.product_indices):
            to_change = []
            for j,reactants_j in enumerate(self.parse.reactant_indices):
                if j not in self.single_molecule_reactions:             # Only include single molecule reactions
                    continue
                for reactanti in reactants_i:
                    if reactanti in reactants_j:
                        to_change.append(j)
            self.dep_graph_I.append(to_change)
            
        # Remove SMRs from dependency graph G (the NRM one)
        for to_change in self.parse.dep_graph:
            to_change_copy = copy.copy(to_change)        
            for reaction in to_change_copy:
                if reaction in self.single_molecule_reactions:
                    to_change.remove(reaction)

    def InitializeTauHeap(self): 
        """ Monte Carlo step to determine all taus and to create an indexed priority queue """           
        # Pre-generate for each reaction random numbers, only used in this function.
        randoms = np.random.random(self.n_reactions)          
        randoms_log_init = -1 * np.log(randoms)      
        self.sim_taus = randoms_log_init/self.sim_a_mu + self.sim_t           # Absolute waiting time for each reaction
        
        # For each reaction find waiting time. For molecular reactions, create list of all reaction times 'on the fly'.
        pairs = []
        self.molecular_taus = {}
        for j in range(self.n_reactions):      
            if j in self.single_molecule_reactions:                
                pair = self.InitializeSingleMoleculeReactions(j)
            else:                                                             # Regular, exponential reaction              
                pair = (j, self.sim_taus[j])                                  #(key, priority) pairs                
            pairs.append( pair )
        self.heap = PQDict(pairs)        
        
        # Pre-generate randoms
        self.randoms_log = -1 * np.log(np.random.random(1000))
        self.count  = 0

    def InitializeSingleMoleculeReactions(self, j):
        """ Initialize single molecule reactions """        
        if self.parse.reaction_orders[j] == 0:    
            tau = self.distr_functions[j](*self.distr_parameters[j])        # Generate just one value, no molecular_taus[j] list.
            pair = (j, tau + self.sim_t)
            
        elif self.parse.reaction_orders[j] == 1:
            reactant_index = self.parse.depends_on[j][0] # reactant_indices
            n_reactants = self.X_matrix[reactant_index]
            if n_reactants == 0:
                self.molecular_taus[j] = []  
                pair = (j, np.inf)                                          # Init heap with an infinite.
            else:
                tau_array = self.distr_functions[j](*self.distr_parameters[j], size = n_reactants)
                
                # Check for negative waiting times.
                assert (tau_array >= 0).all(), "Error: Negative waiting time(s) generated by the distribution '{0:s}'.".format(self.distr_functions[j].__name__)

                tau_array += self.sim_t                                     # Makes time absolute (for sequential reaction, then sim_t!=0)
                tau_array.sort()
                self.molecular_taus[j] = tau_array.tolist()
                pair = (j, self.molecular_taus[j].pop(0))                   # Minimum of the taus goes into the heap, molecular_taus remembers the rest of the taus (for the other molecules).                        
        else:
            print("Error: the fSMM does not support 'single molecule reactions' that are second (or higher) order reactions. Please use the SMM method).")
            sys.exit()
        return pair
            
    def RunExactTimestep(self,settings):
        """ Perform a direct SSA time step and pre-generate M random numbers """
        minimum = self.heap.peek()                                # peek() shows item at top of heap (smallest tau)
        self.reaction_index = minimum[0]                          # Pick reaction to executeO(1)
        self.sim_tau = minimum[1]                                 # Pick tau O(1)        
        if self.sim_tau < settings.endtime:
            self.sim_t = self.sim_tau                             # New time                 
            try:
                self.X_matrix += self.N_matrix_transpose[self.reaction_index]         
            except MemoryError as ex:
                print(ex)
                sys.exit()       
        elif self.sim_tau == np.inf:
            self._IsExhausted = True #=break 18-11-2013            
        else: 
            self.sim_t = settings.endtime
            self.reaction_index = np.nan
            
    def UpdateSystem(self): #19-11-2013
        """ Updates propensities and the tau heap. """           
        ##0 Remove/replace single molecule reaction, if fired.
        if self.reaction_index in self.single_molecule_reactions:  
            if self.parse.reaction_orders[self.reaction_index] == 0: 
                tau_new = self.distr_functions[self.reaction_index](*self.distr_parameters[self.reaction_index]) + self.sim_t
                assert tau_new >= 0, "Error: Negative waiting time(s) generated by the distribution '%s'." % self.distr_functions[self.reaction_index].__name__
                self.heap.updateitem(self.reaction_index, tau_new)
            else:
                # Just remove tau from priority queue
                self.ReplaceMolecularTau(self.reaction_index) 
                
        ##1 Remove single molecule consumed
        ##2 Add Single Molecule Products
        self.UpdateSingleMoleculeReactions()
        
        ##4 Update NRM propensities and waiting times.
        self.UpdateNRM()
        
    def UpdateSingleMoleculeReactions(self):
        """ Remove waitingtimes for molecules reacted, add waitingtimes for molecules formed. """        
        # Remove waiting times
        for j in self.dep_graph_H[self.reaction_index]:               
            if j == self.reaction_index:    # Note: combinable with deleting when own reaction fired?
                continue                    # Already deleted by "if self.reaction_index in self.single_molecule_reactions:".  
            # Reactant of reaction j has reacted in another reaction, remove one at random.         
            random_molecule_index = np.random.randint(low = -1, high = len(self.molecular_taus[j]) )    
            
            # Molecule currently in the heap has reacted in another reaction. Replace this tau by new tau from list.
            if random_molecule_index == -1: 
                self.ReplaceMolecularTau(j)
            else: #Molecule in the tau list has reacted. Remove it from pending list.
                self.molecular_taus[j].pop(random_molecule_index)   
        
        # Add waiting times
        for j in self.dep_graph_I[self.reaction_index]:
            #Assume: first order reaction, reactant in rate equation and one of the reactant formed in a step by other reaction.
            tau_molecular = self.distr_functions[j](*self.distr_parameters[j]) + self.sim_t
            assert tau_molecular >= 0, "Error: Negative waiting time(s) generated by the distribution '{0:s}'.".format(self.distr_functions[j].__name__)
            
            #Check whether this new molecule is the first to fire for this reaction.
            if self.heap[j] == np.inf: 
                self.heap.updateitem(j, tau_molecular) #Replace the infinite by the *tau_molecular*
            elif tau_molecular < self.heap[j]:
                self.molecular_taus[j].insert(0, self.heap[j]) #Put current pending tau back in molecular_taus (at front)
                self.heap.updateitem(j, tau_molecular)
            else:
                self.molecular_taus[j].append(tau_molecular)   
    
    def UpdateNRM(self): 
        """ Step 5 of Gibson and Bruck NRM algorithm. """
        #Set which propensities should be updated
        self.prop_to_update = self.parse.dep_graph[self.reaction_index]   # Propensities to update 
        
        if self.count >= (1000-len(self.prop_to_update)): # TODO: more robust way of pregenerating random numbers. How many randoms max used in a step?
            # Pre-generate random numbers again
            randoms = np.random.random(1000)              
            self.randoms_log = -1 * np.log(randoms)
            self.count = 0
            
        self.Propensities()                                             # Update Propensities  
        
        for j in self.prop_to_update:
            if j == self.reaction_index or self.sim_a_mu_prev[j] == 0:  # Step 5c Gibson Bruck ##Note: due to zero propensity (and an inf tau), it is elaborate to reuse random. Assuming getting new random is faster.
                tau_new = self.randoms_log[self.count]/self.sim_a_mu[j] + self.sim_t
                self.count += 1
            else: #5b, changed propensity
                # if self.sim_a_mu_prev[j] == 0:                       
                    # tau_new = self.randoms_log[self.count]/self.sim_a_mu[j] + self.sim_t
                    # self.count += 1
                # else:
                tau_alpha = self.sim_taus[j]                            # Faster than getting from heap directly (self.heap[j])
                tau_new = self.sim_a_mu_prev[j]/self.sim_a_mu[j]*(tau_alpha - self.sim_t) + self.sim_t              
            #Note, no exception for self.sim_a_mu == 0. Then tau_new automatically becomes infinte (faster than if statement to except this)
        
            self.sim_taus[j] = tau_new #Keep track of the reaction taus, parallel to the heap. Used in getting tau_alpha.
            self.heap.updateitem(j, tau_new) #Updates the tau with index j and resorts the heap.

    def ReplaceMolecularTau(self,j):
        """ Replace tau in the heap by lowest tau in the (pending) list """
        if self.molecular_taus[j]: 
            self.molecular_taus[j].sort()
            tau_new = self.molecular_taus[j].pop(0) #First element contains first time the reaction will occur
            self.heap.updateitem(j, tau_new)
        else: #No pending taus
            self.heap.updateitem(j, np.inf)            
       
        
    def Propensities(self):
        """ Determines the propensities to fire for each reaction at the current time point. At t=0, all the rate equations are compiled. """   
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
        
        assert self.sim_a_mu.min() >= 0, "Error: Negative propensities are found" 
        self.sim_a_mu = abs(self.sim_a_mu)                        # -0 to 0
        self.sim_a_0 = self.sim_a_mu.sum()        
        
