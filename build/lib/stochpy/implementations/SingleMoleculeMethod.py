#! /usr/bin/env python
"""
Single Molecule Method (SMM)
============================

Performs full single molecule stochastic simulation with support for also second-order reactions.

Note: do not use net stoichiometries

Written by M. Moinat and T.R. Maarleveld, Amsterdam, The Netherlands
Adapted from NextReactionmethod.py by T. Maarleveld.
E-mail: tmd200@users.sourceforge.net
Last Change: June 08, 2015
"""

from __future__ import division, print_function, absolute_import

import sys,copy,time,os,re,functools

from .StochPyTools import __species__,StochPySSA_Shared,np

#############################################################

class SingleMoleculeMethod(StochPySSA_Shared):
    """ 
    Performs full single molecule stochastic simulation with single molecule support for second order reactions.
    
    Input:  
     - *File* filename.psc
     - *dir* /home/user/Stochpy/pscmodels/filename.psc    
    """
    def __init__(self,model_file,model_dir,IsQuiet=False):
        self.Parse(model_file,model_dir,IsQuiet=IsQuiet,IsTauleaping = False, IsNRM = False, IsDelayed = False, IsSMM =True)      
    
    def EvaluatePropensities(self):
        """Evaluate propensities to the rate constants, to get the rate parameter for automatically set exponential reactions. """
        species_amount = 10000 #Replace by 1000, divide by 1000        
        for j in self.auto_exponential_reactions:
            propensity = self.parse.propensities[j]
            # Fill in the current 'species_amount' for fixed species.
            for s_id, s_amount in zip(self.fixed_species,self.fixed_species_amount):
                propensity = re.sub('(__species__.{0:s})([^\w]|$)'.format(s_id), str(s_amount) + r'\2', propensity)     
           
            # Replace all other species by <species_amount>
            propensity, n_subs = re.subn('__species__.\w+', str(species_amount), propensity)           
            
            # Evaluates the propensity formula and divide by <species_amount> to correct for each replaced in species.
            k = eval(propensity)/(species_amount**n_subs)            
            try:
                k = float(k)
            except ValueError:
                print("Warning: Defaulting to rate=1 for {0:s}".format(propensity) )
                k = 1.0             
            if k == 0:                                        # k = 0, scale = large (infinite)
                self.distr_parameters[j] = np.inf
            else:
                self.distr_parameters[j] = [ 1.0/k ]          # Convert from rate to scale and make iterable for unpacking            

        
    def Execute(self,settings,IsStatusBar=False):
        """
        Generates T trajectories of the Markov jump process.

        Input:
         - *settings* (class object)   
        """
        if settings.IsSeed:
            np.random.seed(5)   
                
        self.settings = settings       
        self.fixed_species_amount = copy.deepcopy(self.parse.fixed_species_amount)
        self.X_matrix = copy.deepcopy(settings.X_matrix)        
        self.sim_t = copy.copy(settings.starttime)
            
        try:
            self.volume_code = settings.volume_code
        except AttributeError: 
            self.volume_code = "self._current_volume = 1"    # No volume_code present in settings (normal mode)
        exec(self.volume_code)           
             
        self.EvaluatePropensities()         # 31-03-2014. Evaluate propensities of automatically set exponential reactions.
        self.InitializeTauArrays()       
            
        if not self.sim_t:          # Only at sim_t == 0
            self.timestep = 1
            self.sim_output = []
            self.propensities_output = []              
            self.V_output = []        
            self._IsTrackPropensities = False 
            self.SpeciesSelection() 
                 
            self.SetEvents()  
            if not settings.IsOnlyLastTimepoint:
                self.Initial_Conditions()               
        
        nstep_counter = 1
        t1 = time.time()
        while self.sim_t < settings.endtime and self.timestep < settings.timesteps:
            self.n_reactionsPending = functools.reduce(lambda x, y: x + y.size, self.tau_arrays, 0)  # Sums sizes of all arrays of taus.
            if self.n_reactionsPending <= 0:
                settings.endtime = 10**50
                break                                                                       # No molecules can react anymore           
                
            self.RunExactTimestep()            
            self.HandleEvents()
            
            if not self._IsPerformEvent:
                if (self.sim_t < settings.endtime) and (self.timestep < settings.timesteps):     
                    self.UpdateHeap()
                    #print(len(self.tau_arrays[0]),len(self.tau_arrays[1]))
            else:
                # If event performed, recalculate rates (for change in fixed species) and completely rebuild tau arrays.                
                self.EvaluatePropensities()   
                self.InitializeTauArrays()
                
            if not settings.IsOnlyLastTimepoint:
                self.GenerateOutput()   
            
            if (self.sim_t < settings.endtime):         
                self.timestep += 1
                
            t2 = time.time()
            if IsStatusBar and t2-t1> 1:
                t1 = time.time()
                sys.stdout.write('\r')
                sys.stdout.write('simulating {0:s}'.format('.'*nstep_counter) ) 
                sys.stdout.flush() 
                nstep_counter+=1
                if nstep_counter > 10:
                    nstep_counter = 1
                    
        if settings.IsOnlyLastTimepoint or settings.endtime != 10**50:               
            self.GenerateOutput()
            
        if IsStatusBar and t1 and not settings.quiet:  
            sys.stdout.write('\rsimulation done!               \n')         
         
    def RunExactTimestep(self):
        """ Perform a direct SSA time step """ 
        # Get minimum tau
        self.minimums = [(taus.min(),j) for j,taus in enumerate(self.tau_arrays) if taus.size > 0 ] 
        minimum = min(self.minimums)
        self.reaction_index = minimum[1]   
        self.tau = minimum[0]       
        
        if self.tau < self.settings.endtime:
            self.sim_t = self.tau                                      # Tau is absolute
            try:
                self.X_matrix += self.N_matrix_transpose[self.reaction_index]                
            except MemoryError as ex:
                print(ex)
                sys.exit()
        else: 
            self.sim_t = self.settings.endtime
            self.reaction_index = np.nan

    def InitializeTauArrays(self):
        """
        Initializes a list of the arrays which will contain for every reaction the absolute time of reaction.
       
        1) zero-order reaction: an array of one element is created.
        2) first-order reaction: an array with one row and n columns is created.
        3) second-order reaction: an array with n rows and m columns is created. 
        4) molecule number 0, an empty array is created.
        """           
        self.tau_arrays = []
        for j in range(self.n_reactions):
            #reactant_indices = self.parse.reactant_indices[j] 
            depends_on = self.parse.depends_on[j]            
            array_shape = [self.X_matrix[i] for i in depends_on]  # Also works for zero order and no reactants.                  
            #print(0,depends_on)
            #print(1,array_shape)
            try:      
                if self.distr_parameters[j] != np.inf:
                    tau_array = self.distr_functions[j](*self.distr_parameters[j], size = array_shape) # use provided distribution and parameters to generate a tau_array
                else:
                    tau_array = np.array([]) # empty array
            except MemoryError:
                print("Error: Too many putative reaction times have to be generated")                
                sys.exit()            
            #print(2,tau_array)
            # Check for negative waiting times.
            assert (tau_array >= 0).all(), "Error: Negative waiting time(s) generated by the distribution '{0:s}'!".format(self.distr_functions[j].__name__)
            tau_array += self.sim_t  # Makes time absolute (for sequential reaction, then sim_t!=0)                               
            
            if self.parse.reaction_orders[j] == 2 and depends_on[0] == depends_on[1]:       # If e.g. y+y>y2, make y1 reacting with itself infinite.
                np.fill_diagonal(tau_array, np.inf)
            
            self.tau_arrays.append( tau_array )
        
    def UpdateHeap(self):       
        """ Updates the heap by removing reacted molecules and adding formed molecules. Draws waiting time for zero order reaction."""
        self.Remove_Taus_Reactants()
        self.Add_Taus_Products()        
        # If zero order reaction fired, replace its tau.
        if self.parse.reaction_orders[ self.reaction_index ] == 0: 
            j = self.reaction_index       
            new_tau = self.distr_functions[j](*self.distr_parameters[j])
            assert new_tau >= 0, "Error: Negative waiting time(s) generated by the distribution '{0:s}'!".format(self.distr_functions[j].__name__)
            self.tau_arrays[j] = new_tau + self.sim_t      

    def Remove_Taus_Reactants(self):
        """        
        For the molecules reacted in the last reaction, delete the waiting times (taus) for these specific molecules.
        """
        # Search for specific molecules that reacted.
        molecules = [x[0] for x in np.where(self.tau_arrays[self.reaction_index] == self.tau)]  # Specific molecule(s) that reacted. If self.tau occurs more than once, the where takes the first.                
        #species_reacted = self.parse.reactant_indices[self.reaction_index]                     # Species' that reacted in last reaction.       
        species_reacted = self.parse.depends_on[self.reaction_index]                            # Species' that reacted in last reaction.         
        molecules_reacted = list(zip(molecules, species_reacted))                               # =[(mol1,sp1),(mol2,sp2)]; sp1_mol1 reacted with sp2_mol2. e.g. x_3 with y_4.        
        molecules_reacted.sort( reverse=True )                                                  # Prevents indexing error. (If y needs to be deleted twice, always delete the one with higher index first. )          
        # Delete the molecules from every array they occur in.
        for (mol_num, s_index) in molecules_reacted:
            reactions_affected = self.parse.species_depends_on[ s_index ]            
            for j in reactions_affected:                                                       
                if self.distr_parameters[j] != np.inf:                                              # reactions with "inf" have k=0                    
                    if self.parse.reaction_orders[j] == 1: 
                        self.tau_arrays[j] = np.delete(self.tau_arrays[j], mol_num)
                    elif self.parse.reaction_orders[j] == 2:                                        
                        reacts = self.parse.depends_on[j]     
                        axis = reacts.index(s_index)                                                # Specifying molecule on row or column
                        self.tau_arrays[j] = np.delete(self.tau_arrays[j], mol_num, axis)           
                        if reacts[0] == reacts[1]:                                                  # e.g. y + y > z. axis was 0.
                            self.tau_arrays[j] = np.delete(self.tau_arrays[j], mol_num, 1)          # Also delete the molecule on the columns.
                    else:                                                                           # An error should be raised earlier, but just to be safe
                        raise Warning("Error: The SMM does not support third-order reactions. Use the fSMM and model it as a normal reaction.")                    
                        
    def Add_Taus_Products(self):
        """ For the molecules produced in the last reaction, add the absolute tau for every reaction they are reactant in. """
        species_produced = self.products[self.reaction_index] # self.parse.product_indices
        for s_index in species_produced:
            reactions_affected = self.parse.species_depends_on[s_index]
            for j in reactions_affected:
                #reacts = self.parse.reactant_indices[j]
                reacts = self.parse.depends_on[j]
                #print(3,j,self.tau_arrays[j])
                if self.parse.reaction_orders[j] == 1: 
                    new_tau = self.distr_functions[j](*self.distr_parameters[j])
                    assert new_tau >= 0, "Error: Negative waiting time(s) generated by the distribution '{0:s}'!".format(self.distr_functions[j].__name__)
                    self.tau_arrays[j] = np.append( self.tau_arrays[j], new_tau + self.sim_t )                
                elif self.parse.reaction_orders[j] == 2:    # Second-order reaction
                    axis = reacts.index(s_index)          # Check whether produced species on row (=0) or column (=1)
                    temp_arr = self.tau_arrays[j]
                    if axis == 0:
                        new_row = self.distr_functions[j](*self.distr_parameters[j], size = temp_arr.shape[1])  # Taus for every molecule it reacts with
                        assert (new_row >= 0).all(), "Error: Negative waiting time(s) generated by the distribution '{0:s}'!".format(self.distr_functions[j].__name__)
                        self.tau_arrays[j] = np.row_stack( (temp_arr, new_row + self.sim_t) )
                    elif axis == 1:
                        new_col = self.distr_functions[j](*self.distr_parameters[j], size = temp_arr.shape[0])
                        assert (new_col >= 0).all(), "Error: Negative waiting time(s) generated by the distribution '{0:s}'!".format(self.distr_functions[j].__name__)
                        self.tau_arrays[j] = np.column_stack( (temp_arr, new_col + self.sim_t) )
                    
                    if reacts[0] == reacts[1]:                                                                  # axis was 0, row already filled
                        new_col = self.distr_functions[j](*self.distr_parameters[j], size = temp_arr.shape[0]+1)
                        assert (new_col >= 0).all(), "Error: Negative waiting time(s) generated by the distribution '{0:s}'!".format(self.distr_functions[j].__name__)
                        new_col[-1] = np.inf                                                                    # Put infinite on the diagonal
                        self.tau_arrays[j] = np.column_stack( (self.tau_arrays[j], new_col + self.sim_t) )        
                else:                                                                           # An error should be raised earlier, but just to be safe
                    raise Warning("Error: The SMM does not support third-order reactions. Use the fSMM and model it as a normal reaction.")
                #print(4,j,self.tau_arrays[j])
