"""
Extrapolates species copy numbers of single lineage data to extant cell population.

Three steps:
 - Sampling the data at fixed intervals
 - Bin the data
 - Calculating the mass density function form the binned data.
 
Written by T.R. Maarleveld and M. Moinat, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: July 31, 2015
"""

from __future__ import division, print_function, absolute_import
import numpy as np,time,sys,scipy

from  collections import Counter

class ExtrapolateExtant:    
    """
    This class allows for the calculation of population statistics from the simulation of a single lineage        
    """       
    
    def get_interdivision_pdfs(self,n_bins_IDT):
        self._fm_hist, fm_edges = np.histogram(self.data_stochsim_celldivision.interdivision_times, n_bins_IDT, density=True)             
        tau = np.array([(x+y)/2. for (x,y) in zip(fm_edges,fm_edges[1:])])  # (=mean of left and right bin edge)        

        self._fb_hist = 0.5 * np.exp(self.sim_population_growth_rate*tau)*self._fm_hist
        self._fe_hist = self._fb_hist * 2*(1-np.exp(-self.sim_population_growth_rate*tau))            
         
    def get_deterministic_volume(self,n_age_bins): # TODO: different name? remove?: see function in stochpycelldiv.
        """
        Generate deterministic volumes   
        """       
        n_generations = len(self.data_stochsim_celldivision.generation_timesteps)     
        volume_at_birth = np.insert(self.data_stochsim_celldivision.volume_at_birth,0,self.sim_initial_volume) # add initial volume
        n=0   
        g=0
        self._sampled_volumes = np.array([])                
        for age in self._sampled_ages:
            if not n%n_age_bins:
                v0 = volume_at_birth[g]
                g+=1            
                       
            if self.sim_growth_type.lower() == 'exponential':            
                self._sampled_volumes = np.append(self._sampled_volumes,v0*np.exp(self.k*age)) # exponential growth
            elif self.sim_growth_type.lower() == 'linear':
                self._sampled_volumes = np.append(self._sampled_volumes,v0+self.k*age) # linear growth:
            else:
                raise Warning("This growth type '{0}' is not supported".format(self.sim_growth_type))                
            n+=1
            
            
    def sample_per_generation_fixed_idt(self,n_age_bins):
        """
        Sample per generation for a simulation with a fixed interdivision time
        
        Input:
          - *n_age_bins* (integer) number of fixed width samples to make.
        """    
        # assume T = fixed        
        t1 = self.data_stochsim_celldivision.time[0]
        tend = self.data_stochsim_celldivision.time[-1]
        current_steps = 0
        n_generations = len(self.data_stochsim_celldivision.generation_timesteps)           
             
        self._sampled_ages = []        
        self._sampled_IDTs = []
        self._sampled_species = []
        self._sampled_volumes = []
        for n,t2 in enumerate(self.data_stochsim_celldivision.division_times):
            T = t2 - t1            
            ### age ###
            self.dage = T/float(n_age_bins)
            sample_time = 0.5*self.dage
            sample_indices = []
            sample_times = []            
            timesteps_current_generation = self.data_stochsim_celldivision.time[current_steps:current_steps+self.data_stochsim_celldivision.generation_timesteps[n]]
            while sample_time < T:                               
                sample_index = timesteps_current_generation.searchsorted(t1+sample_time,side='right') -1 + current_steps                
                sample_indices.append(sample_index)   
                sample_times.append(sample_time)
                sample_time += self.dage                
                       
            t1 += self.data_stochsim_celldivision.interdivision_times[n]  
            current_steps += self.data_stochsim_celldivision.generation_timesteps[n]
            self._sampled_species.append(self.data_stochsim.species[sample_indices])                 
            self._sampled_ages.append(sample_times)                        
            self._sampled_IDTs.append(np.full(len(sample_times),T)) 
            self._sampled_volumes.append(self.data_stochsim.volume[sample_indices])
        self._sampled_species = np.array(self._sampled_species)
                         
        
    def integrate_fixed_idt(self,n_bins_age):
        """
        get species distribution for a sample of extant cells for a fixed interdivision time
        
        Input:
         - *n_bins_age* (integer)
        """   
        self._extant_species_pn = []
        self._extant_species_pn_sum = []
        
        n_generations = len(self.data_stochsim_celldivision.generation_timesteps)        
        sum_indicator_function = len(self._sampled_ages) # this is only allowed when the interdivision time is fixed
        
        integral_output = []
        for i in range(1,n_bins_age+1):                    
            integral,error =  scipy.integrate.quad(lambda a: self.k*np.exp(self.k*(np.log(2)/self.k-a)),self.dage*(i-1),self.dage*i) # TODO: also with linear growth?
            integral_output.append(integral)

        for species_index in range(len(self.data_stochsim.species_labels)):            
            nmax = 0 
            for age_index in range(n_bins_age):
                nmax_ = max(self._sampled_species[:,age_index][:,species_index])
                if nmax_ > nmax:
                    nmax = nmax_
            
            probilities_at_age = []            
            for age_index in range(n_bins_age):                
                species_at_age = self._sampled_species[:,age_index][:,species_index]                                            
                
                bin_counts = np.bincount(species_at_age)                
                if len(bin_counts) < (nmax+1):
                    bin_counts = np.append(bin_counts, np.zeros(nmax+1 - len(bin_counts) ) )    
            
                probilities_at_age.append(list(bin_counts*integral_output[age_index]/n_generations))
            
            copy_number_probabilities = sum(np.array(probilities_at_age))
            
            p = {}
            for n in range(0,nmax+1):
                p[n] = copy_number_probabilities[n]            
            self._extant_species_pn.append(p)             
            self._extant_species_pn_sum.append(sum(copy_number_probabilities))


    def integrate_fixed_idt_volume(self,n_bins_age,n_bins_volume=20): 
        # TODO: this does not function properly
        n_generations = len(self.data_stochsim_celldivision.generation_timesteps)  
        
        Vmax = max(np.array(self._sampled_volumes).flatten()) 
        Vmin = min(np.array(self._sampled_volumes).flatten()) 
        sum_indicator_function = len(self._sampled_ages)
        
        integral_output = []
        for i in range(1,n_bins_age+1):        
            # works for a fixed IDT
            integral,error =  scipy.integrate.quad(lambda a: self.k*np.exp(self.k*(np.log(2)/self.k-a)),self.dage*(i-1),self.dage*i) 
            integral_output.append(integral)              
        
        Vbins = np.linspace(Vmin, Vmax, n_bins_volume)
       
        volume_prob = []            
        for age_index in range(n_bins_age):                
            volume_at_age = np.array(self._sampled_volumes)[:,age_index][:,0]                  
            
            digitized = np.digitize(volume_at_age, Vbins)
            bin_means = [volume_at_age[digitized == i].mean() for i in range(1, len(Vbins))]
            bin_counts = np.bincount(digitized)                                                       
                    
            if len(bin_counts) < len(Vbins):
                bin_counts = np.append(bin_counts, np.zeros(len(Vbins) - len(bin_counts) ) )
            elif len(bin_counts) > len(Vbins):
                bin_counts = np.delete(bin_counts,-1)

            volume_prob.append(list(bin_counts*integral_output[age_index]/n_generations)) # works only for 1 IDT?        
        
        volume_probabilities = sum(np.array(volume_prob))        
        self._extant_volume_pv = {v:volume_probabilities[i]  for i,v in enumerate(Vbins)}           
        self._extant_volume_pv_sum = sum(volume_probabilities)
        
        
    def sample_per_generation(self,n_age_bins):   
        """
        Sample the ages, interdivision times and species at fixed time intervals in each generation
        reason: no shifts in cell ages.
        
        Example:
        t(age=0.5*dage)=, t(age=1.5*dage)=, t(age=2.5*dage)=2, ..., t(age=0.5*dage)=, t(age=1.5*dage)..
        
        Input:
          - *n_age_bins* (integer) number of fixed width samples to make.
        """    
        t1 = self.data_stochsim_celldivision.time[0]
        tend = self.data_stochsim_celldivision.time[-1]
        current_steps = 0
        n_generations = len(self.data_stochsim_celldivision.generation_timesteps) 
             
        self._sampled_ages = []        
        self._sampled_indices = []
        for n,t2 in enumerate(self.data_stochsim_celldivision.division_times):
            T = t2 - t1            
            ### age ###
            self.dage = T/float(n_age_bins) # different dage for different IDTs      
            sample_time = 0.5*self.dage            
            
            timesteps_current_generation = self.data_stochsim_celldivision.time[current_steps:current_steps+self.data_stochsim_celldivision.generation_timesteps[n]]          
            while sample_time < T:                               
                sample_index = timesteps_current_generation.searchsorted(t1+sample_time,side='right') -1 + current_steps                
                self._sampled_indices.append(sample_index)   
                self._sampled_ages.append(sample_time)               
                sample_time += self.dage
                
            ### IDT ###                
            t1 += self.data_stochsim_celldivision.interdivision_times[n]  
            current_steps += self.data_stochsim_celldivision.generation_timesteps[n]
        
        ### create for every time step the corresponding interdivision time
        full_IDT_array = np.zeros(sum(self.data_stochsim_celldivision.generation_timesteps))
        n=0
        for IDT, n_steps in zip(self._interdivision_times, self.data_stochsim_celldivision.generation_timesteps):
            full_IDT_array[n:n+n_steps] = IDT
            n+=n_steps
        
        self._sampled_IDTs = full_IDT_array[self._sampled_indices]
        self._sampled_species = np.array(self.data_stochsim.species[self._sampled_indices])
        self.get_deterministic_volume(n_age_bins)
        
    
    def calculate_extantCellCN(self, n_bins_age,n_bins_IDT, integration_method= 'riemann' ):
        """
        Calculates the probability mass function of the species copy number in extant cells, using the binned data.
        
        Input:    
          - *n_bins_age* (integer)
          - *n_bins_IDT* (integer)     
          - *integration_method* (string) [default = 'riemann']
        """         
        self._extant_species_pn = [] # Stores for each species the p(n) extant
        self._extant_species_pn_sum = []
        for i in range(len(self.data_stochsim.species_labels)):   
            min_amount = self._sampled_species[:,i].min()
            max_amount = self._sampled_species[:,i].max()
            n_bins_species = max_amount - min_amount + 1
            species_values = np.linspace(min_amount,max_amount,n_bins_species)                
            hist3d, hist_edges = np.histogramdd((self._sampled_species[:,i], self._sampled_IDTs, self._sampled_ages), (n_bins_species, n_bins_IDT, n_bins_age))         
            
            # DO ONCE
            IDT_values = [(x+y)/2. for x,y in zip(hist_edges[1],hist_edges[1][1:])]          
            age_values = [(x+y)/2. for x,y in zip(hist_edges[2],hist_edges[2][1:])]     
 
            if len(IDT_values) > 1:
               dIDT = IDT_values[1] - IDT_values[0]
            else:
               dIDT = 1
               
            if len(age_values) > 1:
               dage = age_values[1] - age_values[0]
            else: 
               dage = 1
            
            # Sum over the species, gives 2D array of sum per age and IDT.
            species_sums = hist3d.sum(axis = 0)
            
            # Initialize p
            p = {n:0 for n in species_values}               
            sum_p_n = 0
            if integration_method.lower() == 'riemann': ### 2D Riemann Sum (less accurate) ####                   
                # For every species amount
                for s_amount_hist, n in zip(hist3d, species_values):
                    index_IDT = 0
                    # For every IDT
                    for age_hist, IDT_value in zip(s_amount_hist, IDT_values):
                        index_age = 0
                        # For every age
                        for bin_count, age_value in list(zip(age_hist, age_values)):
                            normalization = species_sums[index_IDT][index_age]                         
                            if normalization > 0: # Then bin_count also > 0. Prevents dividing 0 by 0
                                p[n] += self._fm_hist[index_IDT] * self.h(age_value,IDT_value) * bin_count/normalization * dage * dIDT
                            index_age += 1        
                        index_IDT += 1                
                    sum_p_n += p[n]            
             
            else: ### 2D trapezoidal rule:  http://mathfaculty.fullerton.edu/mathews/n2003/SimpsonsRule2DMod.html                      
                dx = dIDT
                dy = dage
                
                a = IDT_values[0]
                b = IDT_values[-1]            
                c = age_values[0]
                d = age_values[-1]         
                for s_amount_hist,n in zip(hist3d,species_values):
                    f_a_c = 0
                    f_b_c = 0
                    f_a_d = 0
                    f_b_d = 0                    
                    normalization = species_sums[0][0] 
                    bin_count = s_amount_hist[0][0]
                    if normalization > 0:                      
                        f_a_c = self._fm_hist[0] * self.h(c,a) * bin_count/normalization
                
                    normalization = species_sums[-1][0]     
                    bin_count = s_amount_hist[-1][0]
                    if normalization > 0:   
                        f_b_c = self._fm_hist[-1] * self.h(c,b) * bin_count/normalization
                
                    normalization = species_sums[0][-1]                
                    bin_count = s_amount_hist[0][-1]   
                    if normalization > 0:     
                        f_a_d = self._fm_hist[0] * self.h(d,a) * bin_count/normalization
                    
                    normalization = species_sums[-1][-1]        
                    bin_count = s_amount_hist[-1][-1]
                    if normalization > 0:      
                        f_b_d = self._fm_hist[-1] * self.h(d,b) * bin_count/normalization

                    sum_f_x_c = 0
                    sum_f_x_d = 0
                    sum_f_a_y = 0
                    sum_f_b_y = 0
                    sum_f_x_y = 0                
                    x_i = 1                
                    for x in IDT_values[1:-1]:
                        normalization = species_sums[x_i][0]  # x[i],c
                        bin_count = s_amount_hist[x_i][0]         
                        if normalization > 0:      
                            sum_f_x_c += self._fm_hist[x_i] * self.h(c,x) * bin_count/normalization       
                                         
                        normalization = species_sums[x_i][-1] # x[i],d
                        bin_count = s_amount_hist[x_i][-1] 
                        if normalization > 0:      
                            sum_f_x_d += self._fm_hist[x_i] * self.h(d,x) * bin_count/normalization                        
                        x_i += 1
                    
                    y_i = 1                
                    for y in age_values[1:-1]:
                        normalization = species_sums[0][y_i] # a,y[i]
                        bin_count = s_amount_hist[0][y_i]
                        if normalization > 0:      
                            sum_f_a_y += self._fm_hist[0] * self.h(y,a) * bin_count/normalization                    
                            
                        normalization = species_sums[-1][y_i] # b,y[i]
                        bin_count = s_amount_hist[-1][y_i]
                        if normalization > 0:      
                            sum_f_b_y += self._fm_hist[-1] * self.h(y,b) * bin_count/normalization                        
                        y_i += 1                    
                    
                    x_i = 1                
                    for x in IDT_values[1:-1]:     
                        y_i = 1              
                        for y in age_values[1:-1]:  
                            normalization = species_sums[x_i][y_i]
                            bin_count = s_amount_hist[x_i][y_i]  # x[i],y[i]
                            if normalization > 0:      
                                sum_f_x_y += self._fm_hist[x_i] * self.h(y,x) * bin_count/normalization                    
                            y_i += 1
                        x_i += 1
                    #print(f_a_c,f_b_c,f_a_d,f_b_d)
                    #print(sum_f_x_c,sum_f_x_d,sum_f_a_y,sum_f_b_y,sum_f_x_y)                
                    p[n] = 0.25*dx*dy*(f_a_c + f_b_c + f_a_d + f_b_d + 2*(sum_f_x_c + sum_f_x_d + sum_f_a_y + sum_f_b_y) + 4*sum_f_x_y)                    
                    sum_p_n += p[n]
                #print(sum_p_n)
                
            self._extant_species_pn.append(p) 
            self._extant_species_pn_sum.append(sum_p_n)
            
            
    def calculate_extantCellCN_age(self, age, n_bins_IDT, integration_method= 'riemann' ):
        """
        Calculates the probability mass function of the species copy number in extant cells, using the binned data.        

        Input:    
          - *n_bins_age* (integer)
          - *n_bins_IDT* (integer)                     
          - *integration_method* (string) [default = 'riemann']
        """     
        integration_method= 'Riemann'
        if not integration_method.lower() == 'riemann':
            print("To calculate the extant cell volume distribution we use the riemann sum because no other method is implemented yet")       
                      
        extant_species_age_pn = [] # Stores for each species the p(n) extant
        baby_species_age_pn = [] # Stores for each species the p(n) extant
        for species_index in range(len(self.data_stochsim.species_labels)):
            if age.lower() =='birth':
                Arr_species = self.data_stochsim_celldivision.species_at_birth[:,species_index]            
            elif age.lower() == 'division':
                Arr_species = self.data_stochsim_celldivision.species_at_division[:,species_index]
            else:
                print("Error: age should be 'birth' or 'division', not '{0}'. Birth is selected".format(age))
                Arr_species = self.data_stochsim_celldivision.species_at_birth[:,species_index]                

            min_amount = Arr_species.min()
            max_amount = Arr_species.max()
            n_bins_species = max_amount - min_amount + 1
            species_values = np.linspace(min_amount,max_amount,n_bins_species)

            hist2d, hist_edges = np.histogramdd((Arr_species,self._interdivision_times[:-1]),(n_bins_species,n_bins_IDT))  #species_hist
            
            species_sums = hist2d.sum(axis=0)    
            IDT_values = [(x+y)/2. for x,y in zip(hist_edges[1],hist_edges[1][1:])]  
            if len(IDT_values) > 1:
               dIDT = IDT_values[1] - IDT_values[0]
            else:
               dIDT = 1     
               
            p_e = {n:0 for n in species_values}
            p_b = {n:0 for n in species_values}
            for s_amount_hist, n in zip(hist2d, species_values):
                #print(s_amount_hist,n)
                index_IDT = 0
                # For every IDT
                for bin_count, IDT_value in zip(s_amount_hist, IDT_values):                            
                    normalization = species_sums[index_IDT]
                    if normalization > 0: # Then bin_count also > 0. Prevents dividing 0 by 0                            
                        p_e[n] += self._fe_hist[index_IDT] * bin_count/normalization * dIDT        
                        p_b[n] += self._fb_hist[index_IDT] * bin_count/normalization * dIDT 
                    index_IDT += 1  

            extant_species_age_pn.append(p_e)   
            baby_species_age_pn.append(p_b)
        return(extant_species_age_pn,baby_species_age_pn)            
            
            
    def calculate_extantCellCV(self, n_bins_age,n_bins_IDT, n_bins_volume=20, integration_method= 'riemann' ):
        """
        Calculates the probability density function of the cell volume in extant cells, using the binned data.        

         Input:    
          - *age* (age)
          - *n_bins_IDT* (integer)  
          - *n_bins_volume* [default = 20] (integer)            
          - *integration_method* (string) [default = 'riemann']            
        """
        integration_method= 'Riemann'
        if not integration_method.lower() == 'riemann':
            print("To calculate the extant cell volume distribution we use the riemann sum because no other method is implemented yet")                         

        # For volume make a 3D histogram structure and make a list of the index of every bin to be raised.
        # Result: hist3d[volume, IDT, age]                
        Vmin = self._sampled_volumes.min()
        Vmax = self._sampled_volumes.max()       
        volume_values = np.linspace(Vmin, Vmax, n_bins_volume)          
        hist3d, hist_edges = np.histogramdd((self._sampled_volumes, self._sampled_IDTs, self._sampled_ages), (n_bins_volume, n_bins_IDT, n_bins_age))  
        
        # Age and IDT value of each bin (=mean of left and right bin edge).
        IDT_values = [(x+y)/2. for x,y in zip(hist_edges[1],hist_edges[1][1:])]  
        age_values = [(x+y)/2. for x,y in zip(hist_edges[2],hist_edges[2][1:])]  
        
        if len(IDT_values) > 1:
           dIDT = IDT_values[1] - IDT_values[0]
        else:
           dIDT = 1
           
        if len(age_values) > 1:
           dage = age_values[1] - age_values[0]
        else: 
           dage = 1   
        
        # Sum over the volume, gives 2D array of sum per age and IDT.
        self.volume_sums = hist3d.sum(axis = 0)        
        
        p = {v:0 for v in volume_values}                       
        sum_p_v = 0
        for v_value_hist, v in zip(hist3d, volume_values):
            index_IDT = 0
            # For every IDT
            for age_hist, IDT_value in zip(v_value_hist, IDT_values):
                index_age = 0
                # For every age
                for bin_count, age_value in zip(age_hist, age_values):
                    normalization = self.volume_sums[index_IDT][index_age]                         
                    if normalization > 0: # Then bin_count also > 0. Prevents dividing 0 by 0
                        p[v] += self._fm_hist[index_IDT] * self.h(age_value,IDT_value) * bin_count/normalization * dage * dIDT
                    index_age += 1        
                index_IDT += 1                
            sum_p_v += p[v]            

        self._extant_volume_pv = p
        self._extant_volume_pv_sum = sum_p_v
        
        
    def calculate_extantCellCV_age(self, age, n_bins_IDT,n_bins_volume=20,integration_method= 'riemann' ):
        """
        Calculates the probability density function of the cell volume in extant cells, using the binned data.        

         Input:    
          - *age* (string)
          - *n_bins_IDT* (integer)  
          - *n_bins_volume* [default = 20] (integer)            
          - *integration_method* (string) [default = 'riemann']              
        """
        integration_method= 'Riemann'
        if not integration_method.lower() == 'riemann':
            print("To calculate the extant cell volume distribution we use the riemann sum because no other method is implemented yet")   
        
        if age.lower() =='birth':
            Arr_volume = self.data_stochsim_celldivision.volume_at_birth        
        elif age.lower() == 'division':
            Arr_volume = self.data_stochsim_celldivision.volume_at_division[:-1]
        else:
            print("Error: age should be 'birth' or 'division', not '{0}'. Birth is selected".format(age))
            Arr_volume = self.data_stochsim_celldivision.volume_at_birth
        
        Vmin = Arr_volume.min()
        Vmax = Arr_volume.max()       
        volume_values = np.linspace(Vmin, Vmax, n_bins_volume)                  
        
        hist2d, hist_edges = np.histogramdd((Arr_volume,self._interdivision_times[:-1]),(n_bins_volume,n_bins_IDT))  #volume_hist
        
        volume_sums = hist2d.sum(axis=0)    
        IDT_values = [(x+y)/2. for x,y in zip(hist_edges[1],hist_edges[1][1:])]  
        if len(IDT_values) > 1:
           dIDT = IDT_values[1] - IDT_values[0]
        else:
           dIDT = 1
         
        p_e = {v:0 for v in volume_values}
        p_b = {v:0 for v in volume_values}      
        for v_value_hist, v in zip(hist2d, volume_values):            
            index_IDT = 0
            # For every IDT
            for bin_count, IDT_value in zip(v_value_hist, IDT_values):                            
                normalization = volume_sums[index_IDT]
                if normalization > 0: # Then bin_count also > 0. Prevents dividing 0 by 0                            
                    p_e[v] += self._fe_hist[index_IDT] * bin_count/normalization * dIDT
                    p_b[v] += self._fb_hist[index_IDT] * bin_count/normalization * dIDT                
                index_IDT += 1  
        return(p_e,p_b)  
        
                    
    def h(self,a,tau):
        """
        conditional probability of cell age a given that the interdivision time equals tau for a sample of extant cells
        
        - *a* = age (float)
        - *tau* = IDT (float)
        """
        return self.k * np.exp(self.k*(tau-a))
