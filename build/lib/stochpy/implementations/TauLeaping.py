#! /usr/bin/env python
"""
Optimized Tau-Leaping
=====================
This program performs Optimized Explicit Tau-leaping algorithm, which is an approximate version of the exact Stochastic Simulation Algorithm (SSA). Here, an efficient step size selection procedure for the tau-leaping method [1] is used.

[1] Cao. Y, Gillespie D., Petzold L. (2006), "Efficient step size selection for the tau-leaping method", J.Chem. Phys. 28:124-135

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: June 08, 2015
"""

from __future__ import division, print_function, absolute_import

############################## IMPORTS ###################################

import sys,copy,time,os,random
from .StochPyTools import __species__,StochPySSA_Shared,np

########################### END IMPORTS ##################################

class OTL(StochPySSA_Shared):
    """  
    Input:  
     - *model_file* filename.psc
     - *model_dir* /home/user/Stochpy/pscmodels/filename.psc    
    """
    def __init__(self,model_file,model_dir,IsQuiet=False):
        self.Parse(model_file,model_dir,IsQuiet=IsQuiet,IsTauleaping = True) 
         

    def Execute(self,settings,IsStatusBar=False,epsilon = 0.03,critical_reactions=[]):
        """
        Generates T trajectories of the Markov jump process. 
         - *settings* (python class)
         - *epsilon* [default = 0.03] (float)       
         - *critical_reactions* [default = [] ] (list)
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
        except AttributeError: # No volume_code present in settings
            self.volume_code = "self._current_volume = 1"      
            
        self.Propensities(IsTauleaping=True)
        
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
                self.Initial_Conditions(IsTauleaping = True) 

            # Generate tau's for exact or approximate SSA              
            randoms = np.random.random(1000)
            self.randoms_log = np.log(randoms)*-1
            self.randoms = np.random.random(1000)
            self.count = 0    
                
            # Tau-leap settings 
            IsExhausted = False
            self.sim_tauleaping_steps = 0  
            n_exact_steps = 100                                      # perform direct method for "n_exact_steps"                   
            self.sim_epsilon = epsilon
            self.sim_Nc = 10            
            self.g_vector = np.ones(self.n_species)                  # initialize the g-vector
            self.sim_mu  = np.zeros(self.n_species)                  # Init mu  (32.a)
            self.sim_var = np.zeros(self.n_species)                  # Init var (32.b)  
            self._HAS_GVECTOR_OPTIONS = False 
            self._IsNegativeSpecies = False
            
            self._user_defined_critical_reactions = []
            for r in critical_reactions:
                assert r in self.rate_names, "Error: user-defined critical reaction '{0}' does not exist".format(r)
                self._user_defined_critical_reactions.append(self.rate_names.index(r))              
                   
        nstep_counter = 1 
        t1 = time.time()
        while self.sim_t < settings.endtime and self.timestep < settings.timesteps:            
            if not self._IsNegativeSpecies:    	                    # If there are no negative conc. after a Tau-leap step
                self.GetCriticalReactions()
                self.GetG()                
                if self.sim_a_0 <= 0:                               # All reactants got exhausted                      
                    break               
                self.GetMuVar()
                self.GetTauPrime()
            
            ##################### start Possible Feedback loop ##################       
            self.GetMethod()            
            ##### Optimized Tau-Leaping step #####
            if self._IsOTL:
                self.GetTauPrimePrime()
                self.GetK()
                self.Execute_K_Reactions()
                if not self._IsNegativeSpecies: 
                    self.HandleEvents(IsTauleapingStep=True)
                    self.Propensities(IsTauleaping = True)      # update alle propensities
                        
                    if not settings.IsOnlyLastTimepoint:
                        self.GenerateOutput(IsTauleaping=True)                         
                        
                    self.sim_tauleaping_steps += 1
                    self.timestep += 1                
                elif self._IsNegativeSpecies:                       # Start feedback loop                     
                    self.sim_tau_prime /= 2.0   
            elif self._IsExact:                                     # Direct SSA step
                n=1               
                while (n < n_exact_steps) and (self.sim_t < settings.endtime) and (self.timestep < settings.timesteps):                    
                    if self.sim_a_0 <= 0:                           # All reactants got exhausted
                        IsExhausted = True           
                        break
                        
                    self.RunExactTimestep()
                    self.HandleEvents()                                
                                   
                    if self.sim_t < settings.endtime:  
                        if not self._IsPerformEvent:
                            self.species_to_update = self.parse.reaction_affects[self.reaction_index] # Determine vars to update                
                        else:
                            self.species_to_update = [s for s in range(self.n_species)]         

                    self.Propensities(IsTauleaping=False)     # exact SSA step, so no tauleaping stuff
                    
                    if not settings.IsOnlyLastTimepoint:
                        self.GenerateOutput(IsTauleaping=True)

                    n+=1            
                    
                    t2 = time.time()
                    if IsStatusBar and t2-t1> 1:
                        t1 = time.time()                        
                        sys.stdout.write('\rsimulating {0:s}\r'.format('.'*nstep_counter) ) 
                        sys.stdout.flush()
                        nstep_counter+=1
                        if nstep_counter > 10:
                            nstep_counter = 1
                            sys.stdout.write('\rsimulating {0:s}         '.format('.'*nstep_counter) )
                            sys.stdout.flush()
            #################### End possible feedback loop #################
            t2 = time.time()               
            if IsStatusBar and t2-t1> 1:
                t1 = time.time()                
                sys.stdout.write('\rsimulating {0:s}\r'.format('.'*nstep_counter) ) 
                sys.stdout.flush() 
                nstep_counter+=1
                if nstep_counter > 10:
                    nstep_counter = 1
                    sys.stdout.write('\rsimulating {0:s}         '.format('.'*nstep_counter) )
                    sys.stdout.flush()
                    
            if IsExhausted: 
                settings.endtime = 10**50               
                break
                
        if settings.IsOnlyLastTimepoint or settings.endtime != 10**50:         
            self.GenerateOutput(IsTauleaping=True)       

        if IsStatusBar and t1 and not settings.quiet:  
            sys.stdout.write('\rsimulation done!               \n')
            

    def RunExactTimestep(self):
        """ Perform a direct method SSA time step """
        if self.count == 1000:
            randoms = np.random.random(1000)
            self.randoms_log = np.log(randoms)*-1
            self.randoms = np.random.random(1000)
            self.count = 0                       
        
        self.sim_r2  = self.randoms[self.count]                       # Draw random number 2 [0-1]    
        self.sim_tau = self.randoms_log[self.count]/self.sim_a_0      # reaction time generation   
        self.count+=1
 
        if (self.sim_t + self.sim_tau) < self.settings.endtime:
            self.sim_t += self.sim_tau                                # Time update            

            self.reaction_index = 0
            sum_of_as = self.sim_a_mu[self.reaction_index]
            criteria = self.sim_r2*self.sim_a_0
            while sum_of_as < criteria:                               # Use r2 to determine which reaction will occur
                self.reaction_index += 1    	                      # Index
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
   

    def GetCriticalReactions(self):
        """ Determines the critical reactions (as a boolean vector) """
        if not self.sim_t:            
            self._Nt_nan = copy.copy(self.N_matrix_transpose)     
            self._Nt_nan[self._Nt_nan>=0]= np.nan # look only at the metabolites it's going to exhaust            
      
        self.critical_reactions = []        
        output = self.X_matrix.ravel()/abs(self._Nt_nan)
        minima = np.nanmin(output,axis=1)
        for reaction in minima:
            if reaction < self.sim_Nc:
                self.critical_reactions.append(1)
            else:
                self.critical_reactions.append(0)
        for j in self._user_defined_critical_reactions:
            self.critical_reactions[j] = 1
                
    
    def GetG(self):
        """ Determine the G-vector based on options provided in [1], page 6 """ 
        if not self._HAS_GVECTOR_OPTIONS: # Determine once the Gi options metnioned in [1] page 6
            self.options = np.zeros(self.n_species) 
            for i,hor_i in enumerate(self.parse.species_HORs):
                if hor_i == 2:
                    if self.parse.species_max_influence[i] == 1:
                        self.options[i] = 2
                    elif self.parse.species_max_influence[i] == 2:
                        self.options[i] = 3                    
                elif hor_i == 3:
                    if self.parse.species_max_influence[i] == 1:
                        self.options[i] = 4
                    elif self.parse.species_max_influence[i] == 2:
                        self.options[i] = 5
                    elif self.parse.species_max_influence[i] == 3:
                        self.options[i] = 6                   
            self._HAS_GVECTOR_OPTIONS = True  
        
        for i,option in enumerate(self.options):
            if option == 1: self.g_vector[i] = 1
            elif option == 2: self.g_vector[i] = 2
            elif option == 3: self.g_vector[i] = 2 + 1.0/(self.X_matrix[i]-1)            
            elif option == 4: self.g_vector[i] = 3
            elif option == 5:
                try:
                    self.g_vector[i] = 3 + 1.5/(self.X_matrix[i]-1)
                except:
                    self.g_vector[i] = 3
            elif option == 6:
                try: 
                    self.g_vector[i] = 3 + 1.0/(self.X_matrix[i]-1) + 2.0/(self.X_matrix[i]-2)
                except:
                    self.g_vector[i] = 3            


    def GetMuVar(self):
        """ Calculate the estimates of mu and var for each species (i) """        
        for i,v_i in enumerate(self.parse.N_matrix):
            self.sim_mu[i] = np.dot(v_i,self.sim_a_mu)
            self.sim_var[i]= np.dot(v_i*v_i,self.sim_a_mu)            
            

    def GetTauPrime(self):
        """ Calculate tau' """
        part = np.divide(self.X_matrix.ravel(),self.g_vector)*self.sim_epsilon # eps*x[i]/g[i] for all i
        num = np.array([part,np.ones(len(part))])                  # eps*x[i]/g[i] for all i , 1 for all i
        numerator = num.max(axis=0)                                # max(eps*x[i]/g[i],1) for all i
        abs_mu = np.abs(self.sim_mu)                               # abs(mu) for all i
        bound1 = np.divide(numerator,abs_mu)                       # max(eps*x[i]/g[i],1) / abs(mu[i]) for all i
        numerator_square = np.square(numerator)    	
        bound2 = np.divide(numerator_square,self.sim_var)          # max(eps*x[i]/g[i],1)^2 / abs(mu[i]) for all i
        tau_primes = np.array([bound1,bound2])			
        try:
            self.sim_tau_prime = np.min(tau_primes[~np.isnan(tau_primes)]) # min (bound1,bound2)
        except:
            self.sim_tau_prime = 10**6
            
        
    def GetMethod(self,multiple_of_a_0 = 10.0):                      
        """
        Determines for each time step what to perform: exact of approximate SSA 
        
        Input:
         - *multiple_of_a_0* (float) [default = 10.0] default value based on literature [2] (Cao et et. 2006)
        """                 
        criteria = multiple_of_a_0/self.sim_a_0                               
        if self.sim_tau_prime > criteria and self.sim_tau_prime != 10**6:
            self._IsExact = False
            self._IsOTL = True
        else:
            self._IsExact = True
            self._IsOTL = False
            

    def GetA0c(self):
        """ Calculate the total propensities for all critical reactions """
        self.sim_a_0_c = np.dot(self.critical_reactions,self.sim_a_mu)        
        
    
    def GetTauPrimePrime(self):
        """ Calculate Tau'' """
        if self.count == 1000:                                     # Re-generate random numbers
            randoms = np.random.random(1000)  
            self.randoms_log = np.log(randoms)*-1         
            self.count = 0    
        self.GetA0c()  
        if not self.sim_a_0_c:
            self.sim_tau_prime_prime = 10**6
        else:
            self.sim_tau_prime_prime = self.randoms_log[self.count]/self.sim_a_0_c
            self.count+=1    
              

    def GetK(self):        
        """ Determines the K-vector, which describes the number of firing reactions for each reaction. """
        self.K_vector = np.zeros((self.n_reactions,1),dtype = int)
        if self.sim_tau_prime < self.sim_tau_prime_prime:          # tau' < tau''
            self.sim_tau = self.sim_tau_prime            
            for j,is_critical in enumerate(self.critical_reactions):
                if not is_critical:
                    a_j = self.sim_a_mu[j]
                    Lambda = self.sim_tau * a_j       
                    k_j = np.random.poisson(Lambda)                # draw from a Poisson distribution
                    self.K_vector[j] = [k_j]              
        else:
            self.sim_tau = self.sim_tau_prime_prime                # tau' > tau''            
            probabilities = []
            IsCriticalReactions = False
            for j,is_critical in enumerate(self.critical_reactions):         
                a_j = self.sim_a_mu[j]
                if is_critical:
                    IsCriticalReactions = True
                    p = float(a_j)/self.sim_a_0          
                    probabilities.insert(j,p)
                    if p == 1:                                     # Only one critical reaction
                        self.K_vector[j] = [1]
                elif not is_critical:
                    probabilities.insert(j,0.0)
                    Lambda = self.sim_tau * a_j
                    k_j = np.random.poisson(Lambda)
                    self.K_vector[j] = [k_j]                
            if IsCriticalReactions:                                # Bug fixed jan 15 2011
                (prob,index) = GetSample(probabilities)            # Select one crit.reaction that fires once
                self.K_vector[index] = [1]


    def Execute_K_Reactions(self):
        """ Perform the determined K reactions """
        if (self.sim_t + self.sim_tau) < self.settings.endtime: 
            self._IsNegativeSpecies = False  
            x_temp  = copy.copy(self.X_matrix)    
            x_temp += np.dot(self.parse.N_matrix,self.K_vector).ravel()     
            minimal_amount = min(x_temp)
            if minimal_amount < 0:                                 # Check for negatives after the K reactions 
                self.sim_tau = self.sim_tau/2.0
                self._IsNegativeSpecies = True
            else:            
                self.X_matrix = x_temp                             # Confirm the done K reactions                
                self.sim_t += self.sim_tau   
        else:                       
            self.sim_t = self.settings.endtime            
            self.reaction_index = np.nan 

################### Useful functions #########################

def MinPositiveValue(data):
    """
    This function determines the minimum positive value of a provided data list

    Input:
     - *data*
    Output: 
     - *minimum positive value*
    """   
    return min([elem for elem in data if elem > 0])
    

def GetSample(probabilities):  
    """
    This function extracts a sample from a list of probabilities.
    The 'extraction chance' is identical to the probability of the sample.

    Input:
     - *probabilities*: (list)
    Output: 
     - *sample*
     - *sample index*
    """
    output = []
    min_probability = float(MinPositiveValue(probabilities))
    for p in probabilities:
        for i in range(int(100*p/min_probability)):
            output.append(probabilities.index(p))
    random.sample(output,1)
    index = random.sample(output,1)[0]
    return (probabilities[index],index)
