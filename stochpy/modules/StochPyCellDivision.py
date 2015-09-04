 #! /usr/bin/env python
"""
StochPy Cell Division Module
============================

Example of a sequential simulator. This simulator tracks one cell for N generations. 

Most of the functionalities of the SSA module such as plotting and writing to a file are also available in this module.

Written by T.R. Maarleveld and M. Moinat, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: August 10, 2015
"""
from __future__ import division, print_function, absolute_import

############################## IMPORTS ###################################
import os,sys,copy,time,random

try:
    import pickle
except ImportError:
    import Cpickle as pickle

from . import StochSim

from .StochPyCellDivisionPlot import *
from .PyscesMiniModel import IntegrationStochasticDataObj,RegularGridDataObj

from ..tools.Progress_bar import Progress_bar
from ..tools.ParseDistributions import ParseDistributions,convertInput2IndicesAndValues,MakeDistributionFunction
from ..tools.kSolver import k_solver
from ..tools.ExtrapolateExtant import ExtrapolateExtant


class RegularGridDataObj(RegularGridDataObj):
    HAS_AVERAGE_SPECIES_EXTANT_DISTRIBUTIONS = False    
    
    def setSpeciesExtantDistributionAverage(self,mean,std):
        """
        Set means and standard deviations of species data
        
        Input:
         - *mean* (list)
         - *std* (list)
        """
        self.species_extant_distributions_means = mean
        self.species_extant_distributions_standard_deviations = std
        self.HAS_AVERAGE_SPECIES_EXTANT_DISTRIBUTIONS = True   
        

class CellDivision(ExtrapolateExtant,PlottingFunctions):
    """
    Input options:     
     - *model_file*  [default = CellDivision.psc] (string) 
     - *dir* [default = None] (string)
     - *mode* [default = 'generations'] (string)
     - *end* [default = 3] (float) 
     - *IsInteractive* [default = True] (boolean)
     - *IsQuiet* [default = True] (boolean)
    """
    def __init__(self,method='direct',model_file = 'CellDivision.psc',dir = None,mode = 'generations',end=3,IsInteractive=True,IsQuiet=True):       
        if os.sys.platform != 'win32':
            output_dir = os.path.join(os.path.expanduser('~'),'Stochpy',)
            temp_dir = os.path.join(os.path.expanduser('~'),'Stochpy','temp',)
            if dir == None:
                dir = os.path.join(os.path.expanduser('~'),'Stochpy','pscmodels')
        else:
            output_dir = os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy',)
            temp_dir = os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy','temp',)
            if dir == None:
                dir = os.path.join(os.getenv('HOMEDRIVE')+os.path.sep,'Stochpy','pscmodels')

        self.output_dir = output_dir
        self.model_dir = dir       
        self.temp_dir = temp_dir
        self.StochSim = StochSim.SSA(method,model_file,dir,mode='steps',IsTrackPropensities=False,IsQuiet=IsQuiet)        
        self.StochSim._IsCellDivision=True        
        self.sim_mode = mode.lower()
        self.sim_steps = 10**50
        self.sim_time  = 10**50
        if self.sim_mode == 'generations':
             self.sim_generations = end
        elif self.sim_mode == 'time':
             self.sim_time = end
        elif self.sim_mode == 'steps':
             self.sim_steps = end
        else: 
             self.sim_generations = 3
             
        self._IsBetaDistribution = False
        self._IsAnalyzedExtant = False
        self._IsModeSetByUser = False
        self._IsEndSetByUser = False
        
        if IsInteractive:
          try: 
              Analysis.plt.ion()                        # Set on interactive pylab environment
          except Exception as er:
              print(er)
        else:
          try: 
              Analysis.plt.ioff()                       # Set on interactive pylab environment
          except Exception as er:
              print(er)                     
       
        self.HasCellDivisionParameters = np.zeros(5)    # Collects of each setting whether it has been set.
        self.SetDefaultParameters()
        self._gene_duplications = {}
        

    def Method(self, method):
        """
        Set stochastic simulation algorithm to be used.
        
        Input:
         - *method* (string)
        """
        self.StochSim.Method(method)
        self.data_stochsim = None
        self.data_stochsim_grid = None
        self.data_stochsim_celldivision = None  
              

    def Timesteps(self,s):
        """        
        Set the number of time steps to be generated for each trajectory
        
        Input:
         - *s* (integer)
        """      
        try:
            self.sim_steps = abs(int(s))            
            self.sim_generations = 10**50
            self.sim_time = 10**50
            self.sim_mode = 'steps'
            if not self.StochSim._IsQuiet:
                print("Info: The number of time steps is: {0:d}".format(self.sim_steps) )
            self._IsEndSetByUser = True
            self._IsModeSetByUser = True
        except ValueError:
            raise ValueError("The number of time steps must be an integer")
            

    def Endtime(self,t):
        """        
        Set the end time of the exact realization of the Markov jump process
        
        Input:
         - *t* (float)       
        """    
        try:
            self.sim_time = abs(float(t))
            self.sim_generations = 10**50
            self.sim_steps = 10**50            
            self.sim_mode = 'time' 
            if not self.StochSim._IsQuiet:         
                print("Info: The simulation end time is: {0:f}".format(self.sim_time) )
            self._IsEndSetByUser = True
            self._IsModeSetByUser = True
        except ValueError:
            raise ValueError("The end time must be an integer or a float")      
                 
        
    def Generations(self,g):
        """
        Set the end time of the exact realization of the Markov jump process
        
        Input:
         - *t* (float)  
        """
        try:
            self.sim_generations = int(abs(float(g)))
            self.sim_steps = 10**50
            self.sim_time = 10**50            
            self.sim_mode = 'generations'
            if not self.StochSim._IsQuiet:
                print("Info: The numbers of generations is set to: {0:d}".format(self.sim_generations) )
            self._IsEndSetByUser = True
            self._IsModeSetByUser = True
        except ValueError:
            raise ValueError("The end time must be an integer or a float")
            
    
    def Mode(self,sim_mode='generations'):
        """
        Set mode for a stochastic simulation.
        
        Input:
         - *sim_mode* [default = 'generations'] (string) ['time','steps', 'generations']
        """
        self.sim_mode = sim_mode.lower()        
        if self.sim_mode not in ['generations','steps','time']:
            print("*** WARNING ***: Mode '{0}' not recognized using: 'generations', thus 3 generations will be modeled".format(sim_mode) )
            self.sim_mode = 'generations'
            self.sim_generations = 3
            self.sim_steps = 10**50
            self.sim_time = 10**50        
        self._IsModeSetByUser = True   
                  
    
    def Trajectories(self,n):
        """
        Set the number of trajectories to be generated
              
        Input:
         - *n* (integer)
        """
        self.StochSim.Trajectories(n)
        
    
    def Model(self,model_file,dir = None,set_default=True): 
        """
        Give the model, which is used to do stochastic simulations on

        Input:
         - *model_file* filename.psc (string)
         - *dir* [default = None] the directory where File lives (string)
         - *set_default* [default = True] (boolean)
         
        *set_default* is used to reset the cell division parameters 
        """   
        self.StochSim.Model(model_file,dir)
        self.HasCellDivisionParameters = np.zeros(5)
        if set_default:
            self.SetDefaultParameters()
            if not self.StochSim._IsQuiet:
                print('Info: Please mind that the cell division parameters are reset to the default settings.')
        self._gene_duplications = {}
        
        
    def ChangeParameter(self, parameter, value):
        """
        Change parameter value   
        
        Input:
         - *parameter* (string)
         - *value* (float)
        """
        self.StochSim.ChangeParameter(parameter, value)
        #Reset volume dependencies of propensity formula
        self.SetVolumeDependencies(self._IsVolumeDependent, self._VolumeDependencies, self._SpeciesExtracellular)


    def ChangeInitialSpeciesAmount(self, species, value):
        """
        Change initial species copy number
        
        Input:
         - *species* (string)
         - *value* (float)
        """
        self.StochSim.ChangeInitialSpeciesAmount(species, value)        
        self.SetVolumeDependencies(self._IsVolumeDependent, self._VolumeDependencies, self._SpeciesExtracellular)  


    def ChangeInitialSpeciesCopyNumber(self, species, value):
        """
        Change initial species copy number
        
        Input:
         - *species* (string)
         - *value* (float)
        """
        self.StochSim.ChangeInitialSpeciesCopyNumber(species, value)        
        self.SetVolumeDependencies(self._IsVolumeDependent, self._VolumeDependencies, self._SpeciesExtracellular)  
        
        
    def GetTrajectoryData(self,n=1):
        """ 
        Switch to another trajectory, by default, the last trajectory is directly accessible      
       
        Input:
         - *n* [default = 1] (integer)
        """ 
        self.StochSim.GetTrajectoryData(n) 
        self.data_stochsim = copy.copy(self.StochSim.data_stochsim)
        try:      
            file_in = open(os.path.join(self.temp_dir,'{0}{1:d}_CD.dat'.format(self.StochSim.model_file,n)),'rb')
            self.data_stochsim_celldivision = pickle.load(file_in)
            file_in.close()
        except IOError:
            raise IOError("Trajectory {0:d} does not exist".format(n))  
                                    
    
    def DumpTrajectoryData(self,n):
        """ 
        Input:
         - *n* (integer)
        """ 
        self.StochSim.DumpTrajectoryData(n)        
        if n == self.data_stochsim_celldivision.simulation_trajectory:       
            filename_out = os.path.join(self.temp_dir,'{0}{1:d}_CD.dat'.format(self.StochSim.model_file,n))          
            f = open(filename_out,'wb')            
            pickle.dump(self.data_stochsim_celldivision,f)   
            f.close()
        else:
            print("Error: Trajectory {0} is currently not selected/ does not exist".format(n))
            
            
    def SetSeeding(self,seed=True):
        """        
        Input:
         - *seed* [default=True] (boolean)        
        """
        self.StochSim.SetSeeding(seed)
         

    def SetQuiet(self,quiet=True):
        """
        Input:
         - *quiet* [default=True] (boolean)
        """    
        self._StochSim._IsQuiet = quiet
        
    
    def SetDelayParameters(self, delay_distributions, nonconsuming_reactions = None):
        """       
        Assign the delay input to the SSA.        
       
        Input:
         - *delay_distributions* (dict) with reaction name (or index) as key and distribution as value.
         - *nonconsuming_reactions* (list) [default = None]

        All reactions are assumed to be consuming reactions. Consuming and nonconsuming reactions are defined according to Xiaodong Cai (2007), "Exact stochastic simulation of coupled chemical reactions with delays", J.Phys. Chem. 126:124108.
       
        Example: SetDelayParameters(delay_distributions = {'R1':('fixed',3), 'R2':('gamma',5,1)}, nonconsuming_reactions = ['R2'])
         - Reaction 'R1' will get a delay of 3 seconds and reaction 'R2' a gamma distributed delay with shape=5 and scale=1.
         - Reaction 'R1' will be a consuming reaction, and 'R2' a nonconsuming reaction.

        Value of *delay_distributions* can be any distribution of numpy.random, e.g.:
         - ('gamma', p1,p2) = np.random.gamma(p1,p2) #Where p1 is the shape parameter and p2 the scale parameter.
         - ('exponential', p1) = np.random.exponential(p1) #Where p1 is the scale parameter (NOT the rate).
         - ('uniform', lower, upper) = np.random.uniform(0,lower,upper)        
        """
        self.StochSim.SetDelayParameters(delay_distributions, nonconsuming_reactions)
        
    
    def SetPutativeReactionTimes(self, distributions): 
        """
        Sets the single molecule putative reaction times distribution.
        
        Input:
         - *distributions* (dict) with reaction name (or index) as key and distribution as value.
          
        Value of *distributions* can be any distribution of numpy.random, e.g.:
         - ('gamma', p1,p2) = np.random.gamma(p1,p2) #Where p1 is the shape parameter and p2 the scale parameter.
         - ('exponential', p1) = np.random.exponential(p1) #Where p1 is the scale parameter (NOT the rate).
         - ('uniform', lower, upper) = np.random.uniform(0,lower,upper)
        """   
        self.StochSim.SetPutativeReactionTimes(distributions)
        self._IsManualReactionTime = np.zeros(self.StochSim.SSA.n_reactions)
        for j in list(self.StochSim.putative_reaction_times):
            self._IsManualReactionTime[j] = 1
            
       
    def Reload(self):
        """ Reload the entire model again. Useful if the model file has changed """
        self.StochSim.Reload()
        self.data_stochsim = None
        self.data_stochsim_grid = None
        self.data_stochsim_celldivision = None
        

    def SetDefaultParameters(self):
        """
        Set default parameters: growth rate, growth type, initial volume, and volume distributions        
              
        With default Phi_beta_mean (2) this gives a symmetrical narrow distribution around 2 on [1.5-2.5].
        Beta Distributions make sure that distributions of cell volumes at division and birth do not overlap
        """
        self.SetGrowthFunction(1.0, 'exponential')
        self.SetInitialVolume(1.0)
        self.SetVolumeDistributions( ('beta',5,5), ('beta',2,2), 2)
                       
        self.SetVolumeDependencies(True)        
        self.SetExactDividingSpecies([])
        self.SetNonDividingSpecies([])
        
        
    def SetGrowthFunction(self, growth_rate, growth_type = 'exponential'):
        """
        Define the fixed growth rate and type (exponential and linear growth are supported)
        
        Input:
         - *growth_rate* (float)
         - *growth_type* (string) [default = 'exponential']
        """
        try:
            self.sim_volume_growth_rate = float(growth_rate)
        except ValueError:
            raise ValueError("The growth rate must be an integer or a float")

        if growth_type.lower() == 'exponential':
            self.Vt_codestring = "self._current_volume = {0:f} * np.exp( {1:f}*( self.sim_t - {2:f} ) )" # needs V0, growth rate, t_division
            self.IDT_func = lambda v_start, v_end: np.log(v_end/v_start)/self.sim_volume_growth_rate            
            self.sim_growth_type = 'exponential'
        elif growth_type.lower() == 'linear':
            self.Vt_codestring = "self._current_volume = {0:f} + {1:f} *( self.sim_t - {2:f} )" # needs V0, growth rate, t_division
            self.IDT_func = lambda v_start, v_end: (v_end - v_start)/self.sim_volume_growth_rate
            self.sim_growth_type = 'linear'
        else:
            raise Warning("Growth type '{0}' not recognized, use 'exponential' or 'linear'.".format(growth_type) )
        
        if not self.StochSim._IsQuiet: 
            print( "Info: The {0} growth rate is: {1}".format(self.sim_growth_type, self.sim_volume_growth_rate) )
        self.HasCellDivisionParameters[0] = 1
            
    
    def SetInitialVolume(self, initial_volume):
        """        
        Input:
         - initial_volume (float)         
        """
        try:
            self.sim_initial_volume = abs(float(initial_volume))
            self.HasCellDivisionParameters[1] = 1
        except ValueError:
            raise ValueError("The initial volume must be a positive integer or float")
            
        
    def SetVolumeDistributions(self, Phi, K, Phi_beta_mean = False):
        """
        Set the volume distribution of a sample of mother cells at division (Phi_m) and the volume partitioning ratio (K)
                
        Input:
         - *Phi* (list)
         - *K* (list) 
         - *Phi_beta_mean* (float) only if Phi is a beta distribution, then it is scaled to the random number from *Phi_distr*
         
        Example: SetVolumeDistribution(Phi=('beta',5,5),K=0.5)        
         
        The number drawn from K is multiplied with Phi to get volume of daughter 1: Vbirth = Vdivision * K
        Note that K must be a symmetrical distribution with a mean of 0.5, otherwise one will create a bias for one daughter.
              
        The *Phi* and *K* distributions can be any distribution of numpy.random, e.g.:
         - ('gamma', p1,p2) = np.random.gamma(p1,p2) #Where p1 is the shape parameter and p2 the scale parameter.
         - ('exponential', p1) = np.random.exponential(p1) #Where p1 is the scale parameter (NOT the rate).
         - ('uniform', lower, upper) = np.random.uniform(0,lower,upper)
         Fixed values for *Phi* and *K* are also possible, e.g.:
         - Phi = ('fixed', 2) and K = ('fixed',0.5) creates a lineage with cells that grow from a volume of 1 to 2.
        """
        # Save settings for k solver
        self._Phi_tuple = Phi
        self._K_tuple = K              
        if Phi[0].lower() == 'beta' and K[0].lower() == 'beta':            
            assert Phi_beta_mean != False, "Beta distribution requires setting the *Phi_beta_mean* argument"            
            assert Phi_beta_mean  > 1.0, "The specified *Phi_beta_mean* must be at least 1.0"
            self.Phi_shift = Phi_beta_mean - 0.5
            self.Phi_distribution = lambda x,y: np.random.__dict__['beta'](x,y) + self.Phi_shift
            self.Phi_parameters = Phi[1::] 
            self.K_distribution =  lambda x,y: (np.random.__dict__['beta'](x,y) * (self.Phi_shift-1) + 1)/(self.Phi_shift+1) # When all default on [0.4 - 0.6]
            self.K_parameters = K[1::]            
            self._IsBetaDistribution = True
        else:            
            if Phi[0].lower() == 'beta':
                if self.sim_growth_type != 'exponential' and not self.StochSim._IsQuiet:
                    print("K is not both beta distributed, the specific growth rate will not be calculated.")
                assert Phi_beta_mean != False, "Beta distribution requires setting the *Phi_beta_mean* argument"            
                assert Phi_beta_mean  > 1.0, "The specified *Phi_beta_mean* must be at least 1.0"    
                self.Phi_shift = Phi_beta_mean - 0.5
                self.Phi_distribution = lambda x,y: np.random.__dict__['beta'](x,y) + self.Phi_shift
                self.Phi_parameters = Phi[1::]   
            else:
                if self.sim_growth_type != 'exponential' and not self.StochSim._IsQuiet:
                    print("Phi and K are not both beta distributed, the specific growth rate will not be calculated.")
                self.Phi_distribution, self.Phi_parameters = MakeDistributionFunction(Phi)
            self.K_distribution, self.K_parameters = MakeDistributionFunction(K)
            self._IsBetaDistribution = False           
            
        self.HasCellDivisionParameters[2] = 1
                
 
    def SetVolumeDependencies(self, IsVolumeDependent=True, VolumeDependencies = [], SpeciesExtracellular = []):
        """ 
        Set volume dependency of reactions. The volume updates are done after every reaction firing.
        
        Input:
         - *IsVolumeDependent* (boolean) [default=True]
         - *VolumeDependencies* (list or True)  User specified which reaction are volume dependent and with which order           
         - *SpeciesExtracellular* (list of str) Species that are not inside the cell and should therefore not be divided upon cell division
         
        If *VolumeDependencies* is not specified, the algorithm will try to determine the volume dependency, depending on the order of the reaction (order = number of reactants minus the extracellular species among those reactants)
        If only reaction name is specified, it by default assumes second order reactions (propensity/V)
          e.g.:
            - ['R2',R3'] --> Reactions R2 and R3 are second order reactions
            - [('R2',4), ('R3',3)] --> R2 is fourth order (propensity/V**3) and R3 is third order (propensity/V**2)          
        """        
        assert isinstance(IsVolumeDependent, bool), "IsVolumeDependent argument must be a boolean (True/False)"
        # Save settings for reset of dependencies at .ChangeInitialCopyNumber and .ChangeParameter
        self._IsVolumeDependent    = IsVolumeDependent
        self._VolumeDependencies   = VolumeDependencies
        self._SpeciesExtracellular = SpeciesExtracellular
        
        if IsVolumeDependent:
            if isinstance(SpeciesExtracellular,str): 
                SpeciesExtracellular = [SpeciesExtracellular]
            try:
                self._species_extracellular_indices = [self.StochSim.SSA.species_names.index(x) for x in SpeciesExtracellular]
                self._species_extracellular = [x for x in SpeciesExtracellular]
            except ValueError:
                raise Warning("{0} is not a species name. Try one of: {1}".format(x, self.StochSim.SSA.species_names) )
                
            volume_propensities = self._BuildPropVolumeDependencies(VolumeDependencies, SpeciesExtracellular)
            self.StochSim.SSA.propensities = volume_propensities
            self._IsPropensitiesVolumeDependent = True            
        else:
            self.Reload() # This resets propensities, but it should not remove all other settings!
            if not self.StochSim._IsQuiet:
                print("Info: your model is reloaded and therefore the propensities are reset.")
            self._IsPropensitiesVolumeDependent = False
        self.HasCellDivisionParameters[3] = 1
        
        
    def SetExactDividingSpecies(self, species):
        """
        Set species (as list of names) that are volume dependent, but divide equally (after replication).
        
        Input:
         - *species* (list)
        """
        if isinstance(species,str): 
            species = [species]
        assert isinstance(species,list), "Argument *species* must be a string or a list of strings"
        try:
            self.exact_species_indices = [self.StochSim.SSA.species_names.index(s_id) for s_id in species]
        except ValueError:
            raise Warning("{0} is not a species name. Try one of: {1}".format(s_id, self.StochSim.SSA.species_names) ) 
        
        self.HasCellDivisionParameters[4] = 1
        
    
    def SetNonDividingSpecies(self, species):    
        """
        Set species (as list of names) that are volume dependent, but do not divide.
        
        Input:
         - *species* (list)
        """
        if isinstance(species,str): 
            species = [species]
        assert isinstance(species,list), "Argument *species* must be a string or a list of strings"    
        try:
            self.non_dividing_species_indices = [self.StochSim.SSA.species_names.index(s_id) for s_id in species]
        except ValueError:
            raise Warning("{0} is not a species name. Try one of: {1}".format(s_id, self.StochSim.SSA.species_names) )
            
            
    def SetGeneDuplications(self,species,duplication_times):
        """
        Set species (genes) to duplicate
        
        15 April 2015: duplication times are relative times to the next division          
        
        Input:
         - *species* (list)
         - *duplication_times* (list)         
        """
        self.SetExactDividingSpecies(species)
        if isinstance(species,str):
            species = [species]
        assert isinstance(species,list), "Argument *species* must be a string or a list of strings"
        self._gene_duplications = {}
        for s_id,t in zip(species,duplication_times):                
            s_index = self.StochSim.SSA.species_pos[s_id]
            s_amount = self.StochSim.SSA.X_matrixinit[s_index]
            self._gene_duplications[s_id] = {'n':2,'s_t0':int(s_amount), 't_duplication':t}    # hard coded to a duplication
            
    
    def _BuildPropVolumeDependencies(self, reactions_order, extracellular_species):
        """
        Returns propensity formulas with volume dependency inserted.
        Either using user specified orders or by determining the order from number of reactants.
        
        Input: 
         - *reactions_order* (dict) or (list) of tuples with reaction name and the order. e.g. {'R1':2}. Second-order is assumed if no order is specified.
         - *extracellular_species* (list) of species not affected by cell volume
         
        *** for internal use only *** 
        """
        self._HasVolumeDependencies = np.zeros( self.StochSim.SSA.n_reactions )
        propensities = copy.deepcopy(self.StochSim.SSA.parse.propensities) # Load initial propensities from the parse
        if reactions_order:            
            indices_order = convertInput2IndicesAndValues(reactions_order, self.StochSim.SSA.rate_names, default_value = 2)
            for reaction_index, order in indices_order:
                if order >= 2:                                             # Add volume dependency according to /V**(order-1)
                    propensities[reaction_index] = '({0})/self._current_volume**({1})'.format(propensities[reaction_index], order-1)
                    self._HasVolumeDependencies[reaction_index] = 1        
        else:                                                              # Adds the divide by V according to the number of reactants
            for j,r_id in enumerate(self.StochSim.SSA.rate_names):
                all_reagents = self.StochSim.SSA.parse.Mod.__nDict__[r_id]['AllReagents']
                reactants = [x for x in all_reagents if x[1] < 0]
                n_reactants = sum([abs(sp[1]) for sp in reactants])

                # Remove the number of extracellular ligands in reactants:
                order = n_reactants - len( set.intersection(set(reactants), set(extracellular_species)) ) #Number of overlapping species of reactants and extracellular
                if order >= 2: # No or one reactants, so propensity is not volume dependent (and don't change it)
                    # Add volume dependency according to /V**(order-1)                    
                    propensities[j] = '({0})/self._current_volume**({1})'.format(propensities[j], order-1)
                    self._HasVolumeDependencies[j] = 1
            if not self.StochSim._IsQuiet:
                print("Info: The volume dependency is automatically implemented, but note that extracellular species have to be set specifically")
        return propensities
        

    def DoCellDivisionStochSim(self, mode = False, end = False, method = False, trajectories = False, IsTrackPropensities = False, species_selection = None, rate_selection = None, IsOnlyLastTimepoint = False,quiet=False):
        """ 
        Do stochastic simulations with cell growth and cell division.
        
        Input:
         - *mode* [default = 'generations'] (str) ['generations,'time','steps']
         - *end* (integer) [default = None]
         - *method* (str) [default = False]
         - *trajectories* (integer) [default = False]
         - *IsTrackPropensities* (boolean) [default = False]
         - *rate_selection* [default = None] (list) of names of rates to store. This saves memory space and prevents Memory Errors when propensities propensities are tracked
         - *species_selection* [default = None] (list) of names of species to store. This saves memory space and prevents Memory Errors (occurring at ~15 species).
         - *IsOnlyLastTimepoint* [default = False] (boolean)
         - *quiet* [default = False] (boolean)
        """    
        # Make sure that the species that are duplicated are in the exact dividing species. If not, add them!
        for species in self._gene_duplications:
            s_index = self.StochSim.SSA.species_names.index(species)
            if s_index not in self.exact_species_indices:
                self.exact_species_indices.append(s_index)
        
        if self.StochSim._IsQuiet:
            quiet = True
        if not self.HasCellDivisionParameters.all():
            raise AttributeError("Not all cell division Parameters were set.")
        
        if species_selection and isinstance(species_selection,str):   
            species_selection = [species_selection]
        if species_selection and isinstance(species_selection,list): 
            for s_id in species_selection:
                assert s_id in self.StochSim.SSA.species_names,"Species {0} is not in the model or species selection".format(s_id)

        self.StochSim._IsTrackPropensities = IsTrackPropensities       
        if rate_selection and isinstance(rate_selection,str):   
            rate_selection = [rate_selection]
            self.StochSim._IsTrackPropensities = True
        if rate_selection and isinstance(rate_selection,list): 
            for r_id in rate_selection:
                assert r_id in self.StochSim.SSA.rate_names, "Reaction {0} is not in the model or reaction selection".format(r_id)
            self.StochSim._IsTrackPropensities = True
      
        if mode != False:
            self.Mode(sim_mode = mode)             
            self._IsModeSetByUser = False
        elif mode == False and self.sim_mode != 'generations' and not self._IsModeSetByUser:
            self.Mode('generations')
            
        if end != False:            
            if self.sim_mode == 'generations':
                self.Generations(end)
                self._IsEndSetByUser = False
            elif self.sim_mode == 'steps':
                self.Timesteps(end)
                self._IsEndSetByUser = False 
            elif self.sim_mode == 'time':
                self.Endtime(end)
                self._IsEndSetByUser = False
            self._IsModeSetByUser = False    
        elif end == False and self.sim_generations != 3 and not self._IsEndSetByUser:
            self.Generations(3)
         
        if method != False: 
            self.Method(method)
            self.StochSim._MethodSetBy = "DoStochSim"
        elif method == False and self.StochSim.sim_method_name != "Direct" and self.StochSim._MethodSetBy == "DoStochSim":
            self.Method("Direct")
 
        if trajectories != False: 
            self.StochSim.Trajectories(trajectories)
            self.StochSim._IsTrajectoriesSetByUser = False
        elif trajectories == False and self.StochSim.sim_trajectories != 1 and not self.StochSim._IsTrajectoriesSetByUser:
            self.StochSim.Trajectories(1)           
         
        if self.StochSim._IsTrackPropensities and self.StochSim.sim_method_name in ['FastSingleMoleculeMethod','SingleMoleculeMethod']:
            if not quiet: print("*** WARNING ***: Propensities cannot be tracked with the single molecule method")
            self.StochSim._IsTrackPropensities = False
                   
        self.StochSim._IsFixedIntervalMethod = False
        self.StochSim.HAS_AVERAGE = False
        self._IsAnalyzedExtant = False
         
        self.data_stochsim = IntegrationStochasticDataObj()
        self.data_stochsim_grid = RegularGridDataObj()
        self.data_stochsim_celldivision = IntegrationStochasticCellDivisionObj()
        self.StochSim.DeleteTempfiles()
        self.StochSim.data_stochsim_grid = RegularGridDataObj()
        if not quiet:        
            if self.StochSim.sim_trajectories == 1: 
                print("Info: 1 trajectory is generated")
            else:            
                print("Info: {0:d} trajectories are generated".format(self.StochSim.sim_trajectories) )
                print("Info: Time simulation output of the trajectories is stored at {0:s} in directory: {1:s}".format(self.StochSim.model_file[:-4]+'(traj).dat',self.StochSim.temp_dir) )
         
        self.CalculateSpecificGrowthRate(debug = False)
               
        # Delayed Method
        if self.StochSim._IsDelayedMethod and self.StochSim.HAS_DELAY_PARAMETERS:
            # Pass delay parameters to delayed SSA implementation.
            self.StochSim.SSA.distr_functions        = copy.copy(self.StochSim.delay_distributions)
            self.StochSim.SSA.distr_parameters       = copy.copy(self.StochSim.delay_distr_parameters)
            self.StochSim.SSA.reactions_Consuming    = copy.copy(self.StochSim.delayed_consuming)
            self.StochSim.SSA.reactions_NonConsuming = copy.copy(self.StochSim.delayed_nonconsuming)
        elif self.StochSim._IsDelayedMethod: # No delay parameters set.
            raise AttributeError("No delay parameters have been set for the model '{0:s}'. Use the function .SetDelayParameters().".format(self.StochSim.model_file) ) #7-1-2014 exit if no delay parameters            
        
        # Single Molecule Method  #TODO: check if this works
        if self.StochSim._IsSingleMoleculeMethod: 
            if not self.StochSim.HAS_PUTATIVE_REACTION_TIMES:
                raise Warning("No distributions have been set for the model '{0:s}'. Use the function .SetPutativeReactionTimes().".format(self.StochSim.model_file) )
                 
            if self.StochSim.sim_method_name == "FastSingleMoleculeMethod" and 2 in [self.StochSim.SSA.order[j] for j in list(self.StochSim.putative_reaction_times)]:
                if not quiet: print("Info: Second order is not supported by the fast Single Molecule Method. Switching to the full Single Molecule Method.")
                self.StochSim.Method('SingleMoleculeMethod')
            
            voldep_And_manualfiringtime = self._HasVolumeDependencies * self._IsManualReactionTime
            if voldep_And_manualfiringtime.any():
                raise Warning("Non-exponential putative firing times in the Single Molecule Method cannot be volume dependent ({0}).\nChange the volume dependency with .SetVolumeDependencies().".format([self.StochSim.SSA.rate_names[j] for j in range(self.StochSim.SSA.n_reactions) if voldep_And_manualfiringtime[j] >= 1]) )
                
            # Pass delay parameters to Single Molecule SSA implementation.
            self.StochSim.SSA.distr_functions  = copy.copy(self.StochSim.putative_reaction_times)
            self.StochSim.SSA.distr_parameters = copy.copy(self.StochSim.putative_reaction_times_distr_parameters)                    

            #If Single Molecule Method, set exponential distributions to reactions not specified in StochSim.putative_reaction_times
            if self.StochSim.sim_method_name == 'SingleMoleculeMethod':
                self.StochSim.SSA.auto_exponential_reactions = []                
                for j in range(self.StochSim.SSA.n_reactions):
                    if j not in self.StochSim.SSA.distr_functions:                       # Don't replace already assigned distributions.
                        self.StochSim.SSA.distr_functions[j] = np.random.exponential
                        self.StochSim.SSA.distr_parameters[j] = np.nan                   # To be specified at start of simulation (self.SSA.EvaluatePropensities)
                        self.StochSim.SSA.auto_exponential_reactions.append(j)                                             

        ### Do simulation            
        progressBar = Progress_bar(cycles_total = self.StochSim.sim_trajectories, done_msg = 'time')
        for self.StochSim._current_trajectory in range(1,self.StochSim.sim_trajectories+1):            
            self.StochSim.settings = StochSim.SSASettings(x_matrix = self.StochSim.SSA.X_matrixinit,timesteps = self.sim_steps,starttime = 0,endtime = 0,track_propensities = self.StochSim._IsTrackPropensities,species_selection=species_selection,rate_selection = rate_selection,last_timepoint=IsOnlyLastTimepoint, seed=self.StochSim._IsSeed, quiet=quiet)       
                                
            if self.StochSim._IsSeed:                 # necessary to get correct seeding of each generation
                np.random.seed(5)
          
            self._volume_at_division = []
            self._volume_at_birth = []          
            self._volume_at_birth_notselected = []    # The daughter that is not chosen
            self._species_at_division = []
            self._species_at_birth= []
            self._species_at_birth_not_tracked = []
            self._interdivision_times = []                           
            self._generation_timesteps = []
            self._ages = []                           # For every time step the age of the generation.          
                                    
            # Generate first mother cell volume at division.
            self._next_end_volume = self.Phi_distribution(*self.Phi_parameters)        
            
            assert self._next_end_volume >= self.sim_initial_volume, "The initial volume ({0}) is larger than the drawn mother volume. Please choose a smaller initial volume.".format( self.sim_initial_volume)
            
            #Run generations       
            self.StochSim.SSA.timestep = 1   
            self._current_generation = 0 
            self._next_start_volume = self.sim_initial_volume
            while self._current_generation < self.sim_generations:
                self._current_generation += 1
                #0 The Division volume of current cell is already generated at division.
                self._current_start_volume = self._next_start_volume
                self._current_end_volume = self._next_end_volume
                
                #1.1 Calculate Interdivision Time
                self._current_IDT = self.IDT_func(self._current_start_volume, self._current_end_volume)
                assert self._current_IDT >= 0, "Negative interdivision time. Volume at birth = {0}, volume at division = {1}.".format(self._current_start_volume, self._current_end_volume)
                
                #1.2 Do SSA, end = IDT, mode = 'time'
                self._ExecuteVolumeSSA()            
                self._GetAges()                                                     # Save age data of current generation         
               
                endtime_current_generation = self.StochSim.SSA.sim_t
                endtime_complete_generation = self.StochSim.settings.starttime + self._current_IDT
                              
                if endtime_current_generation >= endtime_complete_generation:       # Only divide if new IDT (age) is reached and not at max generations
                    self._interdivision_times.append( self._current_IDT )           # Save IDT of previous generation
                    self._volume_at_division.append( self._current_end_volume )     # This volume is reached at division
                    
                    if self._current_generation < self.sim_generations:
                        self._Divide()                                              # Divide volumes, species and selects daughter to follow
                        self._volume_at_birth.append( self._next_start_volume )     # Next generation is started with this start volume
                        if not self.StochSim.SSA._IsSpeciesSelection:
                            self.StochSim.SSA.sim_output.append(self._output_after_division)
                        else:
                            self.StochSim.SSA.sim_output.append([self._output_after_division[i] for i in self.StochSim.SSA.sim_output_indices])
                        ### Set settings for new simulation ### 
                        self.StochSim.settings.starttime += self._current_IDT
                        self.StochSim.settings.X_matrix = copy.copy(self.StochSim.SSA.X_matrix)                    
                        self.StochSim.SSA._IsPerformEvent = True                    # Division event
                        # Keep V_output in synchrony with sim_output
                        self.StochSim.SSA.V_output.append( self._next_start_volume )# Division happens instantaneous
                else:   
                    break                                                           # Either end_steps or end_time is reached, no division

            if self.sim_mode != 'generations':
                self._interdivision_times.append( self._current_IDT )               # Add last IDT, of not completed generation
            self.StochSim.FillDataStochsim()
            self.StochSim.data_stochsim.setVolume(self.StochSim.SSA.V_output)     
            
            p_v,v = np.histogram(self.StochSim.data_stochsim.volume,bins=20,density=True)  # TODO hard coded number of bins    
            v = np.array([(x+y)/2. for (x,y) in zip(v,v[1:])])  # (=mean of left and right bin edge)
            
            self.StochSim.data_stochsim.setVolumeDistribution([v,p_v])                    
            self.StochSim.data_stochsim.setSpeciesConcentrations(self._species_extracellular_indices)                       
            
            self._GetDeterministicData()
            self.FillDataStochSimCellDivision()
            if self.StochSim.sim_trajectories == 1 and not quiet:
                print("Number of time steps {0:d}, End time {1:f}, Completed Generations: {2:d}".format(self.StochSim.SSA.timestep, self.StochSim.SSA.sim_t, len(self._volume_at_division)))
            elif self.StochSim.sim_trajectories > 1:
                self.DumpTrajectoryData(self.StochSim._current_trajectory)
            progressBar.update(quiet=quiet)
        
        self.StochSim.sim_trajectories_done = copy.copy(self.StochSim.sim_trajectories)
        self.StochSim._IsSimulationDone = True
        self.data_stochsim = copy.copy(self.StochSim.data_stochsim)
        try: 
            self.StochSim.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.StochSim.SSA.rate_names,self.StochSim.plot.plotnum,quiet=quiet)
        except:
            self.StochSim.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.StochSim.SSA.rate_names,quiet=quiet)
        
        if len(self._generation_timesteps) > 3:
            self.AnalyzeExtantCells() # done with default settings
            
            
    def _GetDeterministicData(self):
        """
        Generate deterministic ages, time, volume after generating a trajectory
        
        hard coded to about 10**6, 5*10**6, or 10**7 data points
        
        *** For internal use only ***        
        """
        completed_generations = len(self._volume_at_division)        
        if not completed_generations:          # no generation finished, 1000 for both for a nice plot
            age_steps = 10**3
            generation_steps = 10**3
        elif completed_generations <= 1000:    # 1 < generations <= 1000
            generation_steps = 10**3
            age_steps = 10**6/completed_generations   
        elif completed_generations <= 10.000:  # 1000 < generations <= 10000
            generation_steps = 10**2
            age_steps = 10**6/completed_generations          
        elif completed_generations <= 100.000: # 10000 < generations <= 100000
            age_steps = 5*10**6/completed_generations
            generation_steps = 10**2
        else:                                  # generations < 100000
            age_steps = 10**7/completed_generations
            generation_steps = 10
        
        max_age = max(self._interdivision_times)
            
        self._time_deterministic = []    
        self._volume_deterministic = []    
        self._ages_deterministic = []

        start_volumes = [self.sim_initial_volume] + self._volume_at_birth
        if completed_generations == len(self._interdivision_times):   
            end_times = self._interdivision_times
        elif completed_generations < len(self._interdivision_times):
            endtime_current_generation = self.StochSim.SSA.sim_t    
            end_times = self._interdivision_times[:-1] + [endtime_current_generation - sum(self._interdivision_times[:-1])]
        else:           
           raise Warning("This should not happen!")

        t1 = 0
        for v1,t2 in zip(start_volumes,end_times):            
            generation_times = np.linspace(t1,t1+t2,generation_steps).tolist()
            t1 += t2    
            self._time_deterministic += generation_times # identical to 10 digits or so
            
            generation_ages = np.linspace(0,t2,age_steps*t2/max_age).tolist()
            self._ages_deterministic += generation_ages            
            
            if self.sim_growth_type == 'exponential': 
                generation_volumes = list(v1*np.exp(self.sim_volume_growth_rate*np.linspace(0,t2,generation_steps)))
            elif self.sim_growth_type == 'linear':
                generation_volumes = list(v1 + self.sim_volume_growth_rate*np.linspace(0,t2,generation_steps))
            else:
                raise Warning("This growth type '{0}' is not supported".format(self.sim_growth_type))
            self._volume_deterministic += generation_volumes
            

    def _GetAges(self):
        """
        Calculates cell age from (total) time output.
        
        *** For internal use only ***        
        """
        marker = len(self._ages)                                                        # Number of time steps already processed saved.
        starttime_generation = sum( self._interdivision_times )                         # Time up to start of this generation.
        simt_generation = self.StochSim.SSA.sim_output[marker:]                         # Output of the current generation
        cell_ages = [ output[0] - starttime_generation for output in simt_generation ]  # Calculate cell ages: For every output, select time and subtract start time of the generation
        self._generation_timesteps.append( len(cell_ages) )
        self._ages.extend(cell_ages)
        

    def _ExecuteVolumeSSA(self):
        """
        Executes the SSA until certain end values (end_time or end_steps). 
               
        *** For internal use only ***
        """
        # Reset the volume code with current parameters        
        self.StochSim.settings.volume_code = self.Vt_codestring.format(self._current_start_volume, self.sim_volume_growth_rate, sum(self._interdivision_times))

        self.StochSim.settings.endtime += self._current_IDT
        if self.StochSim.settings.endtime > self.sim_time:
            self.StochSim.settings.endtime = self.sim_time
        
        # Run a generation or until end_steps or end time is reached # hard coded to the Direct (Volume) SSA        
        if self.StochSim.settings.starttime and self.StochSim._IsTrackPropensities:
            self.StochSim.SSA._IsInitial = True
            self.StochSim.SSA.volume_code = self.StochSim.settings.volume_code
            self.StochSim.SSA.Propensities()    
                  
            output_step = self.StochSim.SSA.sim_a_mu.tolist()
            output_step.insert(0,self.StochSim.SSA.sim_t)
            if self.StochSim.SSA._IsRateSelection:
                output_step = [output_step[j] for j in self.StochSim.SSA.rate_output_indices]
            self.StochSim.SSA.propensities_output.append(output_step)             

        for n,species in enumerate(self._gene_duplications):                                
            self.StochSim.SSA.parse.Mod.__eDict__['duplication_{0}'.format(n)] = {'assignments': {species: str(self._gene_duplications[species]['n']*self._gene_duplications[species]['s_t0'])}, 'delay': 0.0, 'name': 'duplication', 'trigger': ' operator.ge(_TIME_,{0})'.format(self.StochSim.settings.starttime+(self.StochSim.settings.endtime-self.StochSim.settings.starttime)*self._gene_duplications[species]['t_duplication']),'tsymb': None}             
        
        if self._gene_duplications != {}:
            self.StochSim.SSA.parse.Mod.InitialiseEvents(self.StochSim._IsQuiet)
            self.StochSim.SSA.__eDict__ = copy.deepcopy(self.StochSim.SSA.parse.Mod.__eDict__)
            if self.StochSim.settings.starttime:
                self.StochSim.SSA.SetEvents()            
        
        self.StochSim.SSA.Execute(self.StochSim.settings,IsStatusBar=False) 
 
        
    def _Divide(self, n = 0):
        """ 
        Divides the volume, chooses daughter to follow and divides the species according to daughter volume.
        
        *** For internal use only ***
        """
        # 2 Divide volume according to K
        V_daughter1 = self._current_end_volume * self.K_distribution(*self.K_parameters)
        V_daughter2 = self._current_end_volume - V_daughter1
        
        # 3 Choose daughter to follow, based on volume. The larger the volume, the larger the probability to be chosen.
        self._next_end_volume = self.Phi_distribution(*self.Phi_parameters) # Volume of mother cell that one of the daughters will grow to.        
        delta_t = self.IDT_func(V_daughter2, self._next_end_volume) - self.IDT_func(V_daughter1, self._next_end_volume) # Positive = volume 1 larger, Negative = volume 1 smaller        
        P_daughter1 = np.exp(delta_t * self.sim_population_growth_rate)/(1+np.exp(delta_t * self.sim_population_growth_rate)) 
        
        is_daughter1 = np.random.binomial(n = 1, p = P_daughter1)   # = 1 (with prob P_daughter1) or 0       
        if is_daughter1:
            self._next_start_volume = V_daughter1
            self._volume_at_birth_notselected.append( V_daughter2 ) # Save the daughter that is not selected
        else:
            self._next_start_volume = V_daughter2
            self._volume_at_birth_notselected.append( V_daughter1 )
        
        # Do not continue if the daughter volume is larger than the next volume at division.
        if self._next_start_volume > self._next_end_volume:          
            assert n < 20, "Cannot find valid end volume for next generation" # solves potential infinite loop
            print("Rejected volume at birth {0}. The cell division SSA is no longer exact.".format(self._next_start_volume) )
            # Try again to divide. Draws new partitioning K and volume at division .
            self._Divide( n + 1 )                           # n counts number of tries to get a valid volume at division
            self._volume_at_birth_notselected.pop()         # Remove last appended (wrong) value
            return None                                     # Breaks from current _Divide()
        
        if not self.StochSim._IsTauleaping:
            self._species_at_division.append(self.StochSim.SSA.sim_output[-1][1:-1])
        else:
            self._species_at_division.append(self.StochSim.SSA.sim_output[-1][1:]) # we do not store trajectories
        
        # 4.1 Divide species        
        volume_fraction = self._next_start_volume / self._current_end_volume
        # divide each species between two cells, weighted by volume           
        for i,species_amount in enumerate(self.StochSim.SSA.X_matrix):            
            if i in self._species_extracellular_indices:    # species is not in the cell and should not divide
                self.StochSim.SSA.X_matrix[i] = species_amount
            elif i in self.non_dividing_species_indices:    # non-dividing intracellular species
                self.StochSim.SSA.X_matrix[i] = species_amount             
            elif i in self.exact_species_indices:           # exact partitioning                                 
                if species_amount:
                    if species_amount%2:
                       self.StochSim.SSA.X_matrix[i] = np.floor(species_amount/2) + np.random.binomial(1,0.5)   
                    else:
                       self.StochSim.SSA.X_matrix[i] = species_amount/2
            else:
                if species_amount:                          # value is non-zero, to avoid errors
                    self.StochSim.SSA.X_matrix[i] = np.random.binomial(species_amount, volume_fraction)  

        # Replace last time point with species amounts after division in StochSim.SSA
        self._output_after_division = self.StochSim.SSA.X_matrix.tolist()  
        self._output_after_division += [amount for amount in self.StochSim.SSA.fixed_species_amount]            
        if self.StochSim.SSA.__aDict__ != {}:
            self.AssignmentRules()
            self._output_after_division += [value for value in self.StochSim.SSA.assignment_species]
        
        ### species_at_division - species_at_birth ###
        t_division = self.StochSim.settings.endtime
        self._output_after_division.insert(0, t_division )
        if not self.StochSim._IsTauleaping:
            self._output_after_division.append(np.NAN)      # no reaction occurs at cell division        
        
        if self.StochSim.SSA._IsSpeciesSelection:            
            self._species_at_birth.append([self._output_after_division[i] for i in self.StochSim.SSA.sim_output_indices[1:-1]])
        else:
            if not self.StochSim._IsTauleaping:
                self._species_at_birth.append(self._output_after_division[1:-1])
            else:
                self._species_at_birth.append(self._output_after_division[1:]) 
                
        # 4.2 Divide delayed species
        if self.StochSim._IsDelayedMethod:
            pending_delayed = copy.copy(self.StochSim.SSA.Tstruct[1:-1]) # All pending delayed reactions at cell division (excluding the 'inits' at start and end)
            if len(pending_delayed) != 0:              
                number_inherited = np.random.binomial(len(pending_delayed), volume_fraction)                
            else:
                number_inherited = 0
            inherited_delayed = random.sample(pending_delayed, number_inherited)  # Sampling without replacement
            #print('mother delayed:', len(pending_delayed))
            #print('daughter delayed:', len(inherited_delayed))
            self.StochSim.SSA.Tstruct = [(0, np.nan)] + inherited_delayed + [(np.inf, np.nan)]
            self.StochSim.SSA.Tstruct.sort()
            
   
    def CalculateSpecificGrowthRate(self, n_trapezoidal = 20, debug = False):
        """
        Calculates the specific growth rate from the mother volume distribution (Phi) and partitioning distribution (K).
        
        The specific growth rate is retrieved by numerically solving equation 20 of Painter and Marr [1] for the specific growth rate of the population. For a fast calculation, the trapezoidal rule is used for integration and the Secant method for solving the equation.
        
        [1] Painter P.R. and Marr A.G. (1968), "Mathematics of microbial populations", Annu. Rev. Microbiol. 22:519-548.
        
        Input: 
         - *n_trapezoidal* [default = 20] (integer)
         - *debug* (boolean) [default = True] With debug True, the script checks calculated specific growth rate. This takes some time. For a 'silent' calculation, set debug False.
        """          
        # If exponential growth, then specific growth rates equals the volume growth rate. No calculation is needed.         
        if self.sim_growth_type == 'exponential':
            self.sim_population_growth_rate = self.sim_volume_growth_rate
            return self.sim_population_growth_rate
        
        # Other growth types, calculate specific growth rate.
        if self._IsBetaDistribution:                   
            print("Info: Starting calculation of specific growth rate...")
            sys.stdout.flush()
            
            self._Solver = k_solver()
            self._Solver.set_model(distr_Phi = self._Phi_tuple, distr_K = self._K_tuple, shift = self.Phi_shift,
                             growth_rate = self.sim_volume_growth_rate, growth_type = self.sim_growth_type)
            self._Solver.set_integration_param(n_trapezoidal, doquad = False, lim_quad = 8)              
                        
            #Solver.Plot_all_distr()                 
            if debug:
                self.sim_population_growth_rate = self._Solver.Solve_k()
            else:
                t1 = time.time()
                self.sim_population_growth_rate = self._Solver.Get_k(init_guess = self.sim_volume_growth_rate)
                t2 = time.time()
                if not self.StochSim._IsQuiet:
                    print("Time taken to calculate the specific growth rate:", t2-t1)                
        else:
            print("Info: Cannot calculate specific growth rate without beta distributions.")
            print("Info: Setting specific_growth_rate = growth_rate.")
            self.sim_population_growth_rate = self.sim_volume_growth_rate
                                  
 
    def AnalyzeExtantCells(self, n_bins_age = 15, n_bins_IDT = 15, n_bins_volume = 20, integration_method= 'riemann',lb =0.95,ub=1.05):
        """ 
        Obtains the statistics for a population of extant cells for which we use data binning and numerical integration.
        
        The right number of bins is very important for an accurate numerical integration. This function returns a warning message if the sum of probabilities (which should sum to 1) is outside the provided lower and upper bounds.

        Input:
          - *n_bins_age* [default = 15] (integer)
          - *n_bins_IDT* [default = 15] (integer)
          - *n_bins_volume* [default = 20] (integer)
          - *integration_method*  [default = 'riemann'] (string) ['riemann,'trapezoidal']
          - *lb* (float) lower accuracy bound
          - *ub* (float) upper accuracy bound
        """  
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation."        
        assert self.sim_generations > 1, "Simulate for more than 1 trajectory"
        
        if len(np.unique(self._interdivision_times)) == 1 and n_bins_IDT != 1:        
            n_bins_IDT = 1
            if not self.StochSim._IsQuiet:    
                print("Info: number of interdivision time bins is set to 1")        
        
        # Calculate k, if not calculated already.
        try:
            self.k = self.sim_population_growth_rate
        except AttributeError:           
            self.CalculateSpecificGrowthRate()
            self.k = self.sim_population_growth_rate
            
        if not self.StochSim._IsQuiet:
            print("Info: Starting the analysis of the extant cell population...")
            sys.stdout.flush()
        for n in range(1,self.StochSim.sim_trajectories_done+1):
            if self.StochSim.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)           
                
            ##### ExtrapolateExtant #####
            if n_bins_IDT == 1:       # fixed IDT -- samples are all identical, but we have to correct for cell age 
                self.sample_per_generation_fixed_idt(n_bins_age )        
                self.integrate_fixed_idt(n_bins_age)
                
                #self.integrate_fixed_idt_volume(n_bins_age,n_bins_volume,k) # TODO VOLUME                
                self.get_interdivision_pdfs(n_bins_IDT)
                self.sample_per_generation(n_bins_age)                     
                self.calculate_extantCellCV(n_bins_age,n_bins_IDT)
            else:                     # variable IDT -- samples are not identical. We have to correct for cell age and idt
                self.get_interdivision_pdfs(n_bins_IDT)
                self.sample_per_generation(n_bins_age)                                
                self.calculate_extantCellCN(n_bins_age, n_bins_IDT ,integration_method)

                ### Only useful for a variable IDT, otherwise each sample (baby, mother, and extant) is at birth and division identical.
                p_extant,p_baby = self.calculate_extantCellCN_age('birth',n_bins_IDT)
                self.data_stochsim_celldivision.setSpeciesExtantAtBirthDistributions(p_extant,self.StochSim.sim_species_tracked)
                p_extant,p_baby =  self.calculate_extantCellCN_age('division',n_bins_IDT)                
                self.data_stochsim_celldivision.setSpeciesExtantAtDivisionDistributions(p_extant,self.StochSim.sim_species_tracked)
                self.data_stochsim_celldivision.setSpeciesBabyAtDivisionDistributions(p_baby,self.StochSim.sim_species_tracked)
                
                self.calculate_extantCellCV(n_bins_age, n_bins_IDT)
                p_extant,p_baby = self.calculate_extantCellCV_age('birth',n_bins_IDT)
                self.data_stochsim_celldivision.setVolumeExtantAtBirthDistribution(p_extant)
                p_extant,p_baby = self.calculate_extantCellCV_age('division',n_bins_IDT)
                self.data_stochsim_celldivision.setVolumeExtantAtDivisionDistribution(p_extant)
                self.data_stochsim_celldivision.setVolumeBabyAtDivisionDistribution(p_baby)
            
            # Species     
            self.data_stochsim_celldivision.setSpeciesExtantDistributions(self._extant_species_pn,self.StochSim.sim_species_tracked)   
            self.data_stochsim_celldivision.setVolumeExtantDistribution(self._extant_volume_pv)                      
            
            if self.StochSim.sim_trajectories_done > 1:
                self.DumpTrajectoryData(n)   
        if not self.StochSim._IsQuiet:    
            print("Info: Successfully analyzed the extant cell species distribution")
        self._IsAnalyzedExtant = True
        for n in range(1,self.StochSim.sim_trajectories_done+1):
            if self.StochSim.sim_trajectories_done > 1:
                self.GetTrajectoryData(n) 
                
            # For every species: give a warning 
            for i,extant_species_pn_sum in enumerate(self._extant_species_pn_sum):
                if extant_species_pn_sum < lb or extant_species_pn_sum > ub:
                    print("*** WARNING ***: We found a sum of species probabilities {0:0.3f} in the extant cells which is not between the desired bounds ({1},{2}). Our advice is to specify a different number of bins for both age and interdivision times.".format(extant_species_pn_sum,lb,ub)) 
                    break   
                    
    # TODO: remove? rename? whatever                
    def GetAverageSpeciesExtantDistributions(self):
        """ Get average species distributions """      
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"
        assert self._IsAnalyzedExtant, "First analyze the extant cell population"
        assert not self.StochSim._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"
        
        D_distributions = {}
        for s_id in self.StochSim.sim_species_tracked:
            D_distributions[s_id] = {}
        L_distributions_means = []
        L_distributions_standard_deviations = []
        for n in range(1,self.StochSim.sim_trajectories_done+1): 
            if self.StochSim.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            for i,s_id in enumerate(self.StochSim.sim_species_tracked):                
                for m,s_amount in enumerate(self.data_stochsim_celldivision.species_extant_distributions[i][0]):
                    if not s_amount in list(D_distributions[s_id]):
                        D_distributions[s_id][s_amount] = []
                    D_distributions[s_id][s_amount].append(self.data_stochsim_celldivision.species_extant_distributions[i][1][m])
            
        for s_id in self.StochSim.sim_species_tracked:
            L_amount = list(D_distributions[s_id])  # for a given species 
            L_means = []
            L_stds = []   
            for s_amount in L_amount:
                while len(D_distributions[s_id][s_amount]) < (n-1):
                    D_distributions[s_id][s_amount].append(0)
                L_means.append(np.mean(D_distributions[s_id][s_amount]))
                L_stds.append(np.std(D_distributions[s_id][s_amount]))    
            L_distributions_means.append([L_amount,L_means]) 
            L_distributions_standard_deviations.append([L_amount,L_stds])
        self.data_stochsim_grid.setSpeciesExtantDistributionAverage(L_distributions_means,L_distributions_standard_deviations)                      
        

    def GetRegularGrid(self,n_samples=250):
        """
        The Gillespie method generates data at irregular time points. However, it is possible to output data on a fixed regular time grid where the user can specify the resolution of the grid (n_samples).  
        
        We use a higher number of default samples, because we expect more stochasticity due to explicit modeling of cell divisions
       
        Input:
         - *n_samples* [default = 250] (integer)
        """    
        self.StochSim.GetRegularGrid(n_samples)        
        self.StochSim.data_stochsim_grid.setSpeciesConcentrations(self._species_extracellular_indices)
        (self.StochSim.data_stochsim_grid.species_concentrations_means,self.StochSim.data_stochsim_grid.species_concentrations_standard_deviations) = Analysis.GetAverageResults(self.StochSim.data_stochsim_grid.species_concentrations)
        self.data_stochsim_grid = copy.copy(self.StochSim.data_stochsim_grid)           
        
      
    def PrintSpeciesMeans(self,sample='extant'):
        """ 
        Print the means (3 decimals) of each species for the selected trajectory
        
        Input:
         - *sample* (string) [default = 'extant']
        """
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"      
        assert not self.StochSim._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"
        if sample == 'mother':
            means = self.data_stochsim.species_means
        else:            
            if sample == 'baby': 
                print("Not yet implemented for a sample of baby cells")
            assert self._IsAnalyzedExtant, "First analyze the extant cell population"
            means = self.data_stochsim_celldivision.species_extant_means
            
        for s_id in self.data_stochsim.species_labels:   
            mu = means[s_id]         
            if not mu < 0.001:
                print("{0:s}\t{1:0.3f}".format(s_id,mu))   
            else:
                print("{0:s}\t{1:0.3e}".format(s_id,mu))               
        
        
    def PrintSpeciesStandardDeviations(self,sample='extant'):
        """
        Print the means (3 decimals) of each species for the selected trajectory

        Input:
         - *sample* (string) [default = 'extant']        
        """
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"    
        assert not self.StochSim._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"
        if sample == 'mother':
            standard_deviations = self.data_stochsim.species_standard_deviations
        else:
            if sample == 'baby': 
                print("Not yet implemented for a sample of baby cells")      
       
            assert self._IsAnalyzedExtant, "First analyze the extant cell population"            
            standard_deviations = self.data_stochsim_celldivision.species_extant_standard_deviations
        
        print("Species\tMean")
        for s_id in self.data_stochsim.species_labels:   
            sigma = standard_deviations[s_id]         
            if not sigma < 0.001:
                print("{0:s}\t{1:0.3f}".format(s_id,sigma))   
            else:
                print("{0:s}\t{1:0.3e}".format(s_id,sigma))  
                

    def PrintSpeciesDistributions(self,sample='extant'): 
        """ Print obtained species distributions for the selected trajectory """      
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"      
        assert not self.StochSim._IsOnlyLastTimepoint, "Determining statistics is disabled when saving only the last time point"
        if sample == 'mother':
            distributions = self.data_stochsim.species_distributions
        else:
            if sample == 'baby': 
                print("Not yet implemented for a sample of baby cells")      
        
            assert self._IsAnalyzedExtant, "First analyze the extant cell population"            
            distributions = self.data_stochsim_celldivision.species_extant_distributions        
        
        for i,species_dist in enumerate(distributions):
            print("Copy number ({0:s})\tPMF".format(self.StochSim.sim_species_tracked[i]) )
            for m in range(len(species_dist[0])):
                x = species_dist[0][m]
                p_x = species_dist[1][m]
                if not p_x < 0.001:
                    print("{0:d}\t{1:0.3f}".format(x,p_x))
                else:
                    print("{0:d}\t{1:0.3e}".format(x,p_x))


    def Export2File(self,analysis='timeseries',datatype='species', IsAverage = False, directory=None): # TODO
        """
        Write data to a text document     
    
        Input:
         - *analysis* [default = 'timeseries'] (string) options: timeseries, distribution, mean, std, autocorrelation, autocovariance
         - *datatype*  [default = 'species'] (string) options: species, propensities, waitingtimes
         - *IsAverage* [default = False] (boolean)   
         - *directory* [default = None] (string)
        """        
        self.StochSim.Export2File(analysis,datatype,IsAverage,directory)
        
        
    def ExportExtant2File(self,analysis='distribution',IsAverage = False, directory=None):
        """
        Write data to a text document     
    
        Input:
         - *analysis* [default = 'distribution'] (string) options: distribution, mean, std, autocorrelation, autocovariance        
         - *IsAverage* [default = False] (boolean)   
         - *directory* [default = None] (string)
        """        
        assert self.StochSim._IsSimulationDone, "First do a stochastic simulation"   
        assert not self.StochSim._IsOnlyLastTimepoint, "Plotting is disabled when saving only the last time point"
        assert self._IsAnalyzedExtant , "First analyze the extant cell population"        
        
        if directory == None:
            if not IsAverage:
                directory = os.path.join(self.output_dir,"{0:s}_{1:s}".format(self.StochSim.model_file,analysis))
            else:
                directory = os.path.join(self.output_dir,"{0:s}_{1:s}_{2:s}".format(self.StochSim.model_file,"average_extant",analysis))
        else:
            if not os.path.exists(directory):
                os.makedirs(directory)
            if not IsAverage:
                directory = os.path.join(directory,"{0:s}_{1:s}".format(self.StochSim.model_file,analysis))  
            else:
                directory = os.path.join(directory,"{0:s}_{1:s}_{2:s}".format(self.StochSim.model_file,"average_extant",analysis))
        
        if analysis.lower() in ['distributions','distribution'] and not IsAverage:            
            for n in range(1,self.StochSim.sim_trajectories_done+1):       
                if self.StochSim.sim_trajectories_done > 1:  
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')  
                for i,L_species_dist in enumerate(self.data_stochsim_celldivision.species_extant_distributions):
                    file_out.write("Copy number ({0:s})\tPMF\n".format(self.StochSim.sim_species_tracked[i]) )
                    for m in range(len(L_species_dist[0])):
                        x = L_species_dist[0][m]
                        p_x = L_species_dist[1][m]                        
                        file_out.write("{0:d}\t{1}\n".format(x,p_x))                                   
                file_out.close()
        elif analysis.lower() == 'mean' and not IsAverage: 
            for n in range(1,self.StochSim.sim_trajectories_done+1):       
                if self.StochSim.sim_trajectories_done > 1:  
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')  
                file_out.write("Species\tMean\n") 
                for s_id in self.StochSim.sim_species_tracked: 
                    file_out.write("{0:s}\t{1:f}\n".format(s_id,self.data_stochsim_celldivision.species_extant_means[s_id])) 
                file_out.close()
                if not self.StochSim._IsQuiet:
                    print("Info: Species extant means output is successfully saved at: {0:s}".format(file_path) )
                
        elif analysis.lower() == 'std' and not IsAverage: 
            for n in range(1,self.StochSim.sim_trajectories_done+1):       
                if self.StochSim.sim_trajectories_done > 1:  
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')  
                file_out.write("Species\tStandard Deviation\n")                
                for s_id in self.StochSim.sim_species_tracked: 
                    file_out.write("{0:s}\t{1:f}\n".format(s_id,self.data_stochsim_celldivision.species_extant_standard_deviations[s_id])) 
                file_out.close()
                if not self.StochSim._IsQuiet:
                    print("Info: Species extant standard deviations output is successfully saved at: {0:s}".format(file_path) )           

        elif analysis.lower() in ['distributions','distribution'] and IsAverage:            
            if not self.data_stochsim_grid.HAS_AVERAGE_SPECIES_EXTANT_DISTRIBUTIONS:
                self.GetAverageSpeciesExtantDistributions()
        
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')     
            for i,s_id in enumerate(self.StochSim.sim_species_tracked):
                file_out.write("Copy number\t{0:s} (Mean)\t{0:s} (STD)\n".format(s_id) )   
                for m in range(len(self.data_stochsim_grid.species_extant_distributions_means[i][0])):
                    s_amount = self.data_stochsim_grid.species_extant_distributions_means[i][0][m] 
                    s_probability_mean = self.data_stochsim_grid.species_extant_distributions_means[i][1][m]
                    s_probability_std = self.data_stochsim_grid.species_extant_distributions_standard_deviations[i][1][m]
                    file_out.write("{0:0.0f}\t{1:f}\t{2:f}\n".format(s_amount,s_probability_mean,s_probability_std) )
                file_out.write("\n")    
            if not self.StochSim._IsQuiet:              
                print("Info: Averaged species distributions output is successfully saved at: {0:s}".format(file_path) ) 
        else:
            raise UserWarning("No valid option specified. Nothing is exported. See help function (help(ExportExtant2File))")
            

    def FillDataStochSimCellDivision(self):
        """ Fill the data_stochsim_celldivision object after generating a realization """         
        self.data_stochsim_celldivision.setSimulationInfo(self.StochSim.SSA.timestep,self.StochSim.SSA.sim_t,self.StochSim._current_trajectory,self._generation_timesteps)
        self.data_stochsim_celldivision.setTime(self.StochSim.data_stochsim.time)
        self.data_stochsim_celldivision.setAges(self._ages)
        self.data_stochsim_celldivision.setDeterministicData(self._time_deterministic,self._volume_deterministic,self._ages_deterministic)
        if self._current_generation > 1:
            self.data_stochsim_celldivision.setInterdivisionTimes(self._interdivision_times,self._current_generation)    
            self.data_stochsim_celldivision.setVolumeAtDivision(self._volume_at_division)
            self.data_stochsim_celldivision.setVolumeAtBirth(self._volume_at_birth, self._volume_at_birth_notselected)
            self.data_stochsim_celldivision.setSpeciesAtDivision(self._species_at_division,self.StochSim.sim_species_tracked)
            self.data_stochsim_celldivision.setSpeciesAtBirth(self._species_at_birth,self.StochSim.sim_species_tracked)    
        
            (L_probability_mass, D_means, D_stds,D_moments) = Analysis.GetSpeciesDistributions(self.data_stochsim_celldivision.getSpeciesAtDivision(),self.StochSim.sim_species_tracked)
            self.data_stochsim_celldivision.setSpeciesAtDivisionDistributions(L_probability_mass,D_means,D_stds,D_moments)
            (L_probability_mass, D_means, D_stds,D_moments) = Analysis.GetSpeciesDistributions(self.data_stochsim_celldivision.getSpeciesAtBirth(),self.StochSim.sim_species_tracked)
            self.data_stochsim_celldivision.setSpeciesAtBirthDistributions(L_probability_mass,D_means,D_stds,D_moments)            

            (L_probability_mass, D_means, D_stds,D_moments) = Analysis.GetDataDistributions(self.data_stochsim_celldivision.getVolumeAtDivision(),['v'])
            self.data_stochsim_celldivision.setVolumeAtDivisionDistribution(L_probability_mass,D_means,D_stds,D_moments)
            (L_probability_mass, D_means, D_stds,D_moments) = Analysis.GetDataDistributions(self.data_stochsim_celldivision.getVolumeAtBirth(),['v'])
            self.data_stochsim_celldivision.setVolumeAtBirthDistribution(L_probability_mass,D_means,D_stds,D_moments)   


    def GetWaitingtimes(self):
        """ Get for each reaction the waiting times """ 
        self.StochSim.GetWaitingtimes()
        self.data_stochsim = copy.copy(self.StochSim.data_stochsim)  

            
    def SumExtantSpeciesDistributions(self):
        """ 
        Use this function to test the sum of p(Nx=n). This sum should be (very) close to zero. If not, do AnalyzeExtantCells with a different number of age and IDT bins. 
        """
        for n in range(1,self.StochSim.sim_trajectories_done+1):
            if self.StochSim.sim_trajectories_done > 1:
                self.GetTrajectoryData(n) 
            # For every species
            for i,extant_species_pn_sum in enumerate(self._extant_species_pn_sum):
                print(self.StochSim.sim_species_tracked[i],"{0:0.3f}".format(extant_species_pn_sum))
        
    def ShowSpecies(self):
        """ Print the species of the model """
        self.StochSim.ShowSpecies()
        

    def ShowOverview(self):
        """ Print an overview of the current settings """
        self.StochSim.ShowOverview()            

     
class IntegrationStochasticCellDivisionObj(object):
    """
    This class is specifically designed to store the results of a stochastic time simulation
    It has methods for setting e.g. the Time, Labels, Species, and Volume data and
    getting e.g.  Time, Species,and Volume (including time) arrays. However, of more use:

    - getDataAtTime(time) the data generated at time point "time".
    - getDataInTimeInterval(time, bounds=None) more intelligent version of the above
      returns an array of all data points where: time-bounds <= time <= time+bounds
    """    
    time = None
    species_extant_distributions = None
    species_at_division = None
    species_at_birth = None
    volume_at_division = None
    volume_at_birth = None
    interdivision_times = None
    division_times = None
    ages = None
    generation_timesteps = None    
    species_labels = None
    xdata_labels = None
    time_label = 'Time'
    
    HAS_DETERMINISTIC = False
    HAS_TIME = False    
    HAS_SPECIES_AT_DIVISION = False
    HAS_SPECIES_AT_BIRTH = False
    HAS_AGES = False
    HAS_IDT = False
    HAS_SPECIES_EXTANT = False
    HAS_VOLUME_EXTANT = False
    
    HAS_VOLUME_AT_DIVISION = False
    HAS_VOLUME_AT_BIRTH = False    
    HAS_XDATA = False

    IS_VALID = True    
    TYPE_INFO = 'Stochastic'  
    
    def getSpeciesDistributions(self,distributions,species):
        """
        get species distributions
        
        Input:
         - *distributions* (dict)
         - *species* (list)
        """        
        D_means = {}
        D_stds = {}
        D_moments = {}
        L_probability_mass = []
        for i,s_id in enumerate(species):            
            x = np.array(sorted(distributions[i]),dtype=int)                    
            p_x = np.array([distributions[i][x_i] for x_i in x])
            
            mu = (x*p_x).sum()
            mu_sq = (x**2*p_x).sum()
            var =  mu_sq - mu**2
            std = var**0.5
            L_probability_mass.append([x,p_x])     
             
            D_means[s_id] = mu
            D_stds[s_id] = std
              
            D_moments[s_id] = {}
            D_moments[s_id]['1'] = mu
            D_moments[s_id]['2'] = mu_sq
            D_moments[s_id]['3'] = (x**3*p_x).sum()
            D_moments[s_id]['4'] = (x**4*p_x).sum()  
        return L_probability_mass,D_means,D_stds,D_moments


    def getVolumeDistribution(self,distribution):   
        """
        get volume distribution
        
        Input:
         - *distributions* (dict)
         - *species* (list)
        """                                     
        v = np.array(sorted(distribution))
        p_v = np.array([distribution[v_i] for v_i in v])
        
        mu = (v*p_v).sum()
        mu_sq = (v**2*p_v).sum()
        var =  mu_sq - mu**2
        std = var**0.5
          
        D_moments = {}
        D_moments['1'] = mu
        D_moments['2'] = mu_sq
        D_moments['3'] = (v**3*p_v).sum()
        D_moments['4'] = (v**4*p_v).sum()          
        return [v,p_v],mu,std,D_moments   
    
    
    def setLabels(self, species):
        """
        Set the species labels
        
        Input:
         - *species* a list of species labels
        """
        self.species_labels = species
        

    def setTime(self, time, lbl=None):
        """
        Set the time vector

        Input:
         - *time* a 1d array of time points
         - *lbl* [default=None] is "Time" set as required
        """
        self.time = time #.reshape(len(time), 1)
        self.HAS_TIME = True
        if lbl != None:
            self.time_label = lbl
            

    def setSimulationInfo(self,timesteps,endtime,simulation_trajectory,steps_per_generation):
        """
        set Simulation Information
        
        Input:
         - *timesteps* (integer)
         - *endtime* (float)
         - *simulation_trajectory* (integer)
        """
        self.simulation_timesteps = timesteps
        self.simulation_endtime = endtime
        self.simulation_trajectory = simulation_trajectory
        self.generation_timesteps = np.array(steps_per_generation)
                
    
    def setInterdivisionTimes(self,idt,n_generations):
        """
        set the (inter)division times array
        
        Input:
         - *idt* (list)
        """
        self.interdivision_times = np.array(idt)
        self.division_times = np.cumsum(self.interdivision_times , dtype = 'd')        
        HAS_IDT = True

        
    def setAges(self,ages):
        """
        set the ... ages array
        
        Input:
         - *ages* (list)
        """
        self.ages = np.array(ages)
        HAS_AGES = True
        
        
    def setSpeciesExtantDistributions(self,distributions,species):
        """
        Set the species array
        
        Input:
         - *distributions* (list)
         - *species* (list)
        """
        distributions,means,stds,moments = self.getSpeciesDistributions(distributions,species)
        self.species_extant_distributions = distributions
        self.species_extant_means = means
        self.species_extant_standard_deviations = stds
        self.species_extant_moments = moments
        self.HAS_SPECIES_EXTANT = True 
         
        
    def setSpeciesExtantAtBirthDistributions(self,distributions,species):
        """
        Set the species array
        
        Input:
         - *distributions* (list)
         - *species* (list)
        """
        distributions,means,stds,moments = self.getSpeciesDistributions(distributions,species)
        self.species_extant_at_birth_distributions = distributions
        self.species_extant_at_birth_means = means
        self.species_extant_at_birth_standard_deviations = stds
        self.species_extant_at_birth_moments = moments
      

    def setSpeciesExtantAtDivisionDistributions(self,distributions,species):
        """
        Set the species array
        
        Input:
         - *distributions* (list)
         - *species* (list)
        """
        distributions,means,stds,moments = self.getSpeciesDistributions(distributions,species)
        self.species_extant_at_division_distributions = distributions
        self.species_extant_at_division_means = means
        self.species_extant_at_division_standard_deviations = stds
        self.species_extant_at_division_moments = moments 
        
        
    def setSpeciesBabyAtDivisionDistributions(self,distributions,species):
        """
        Set the species array
        
        Input:
         - *distributions* (list)
         - *species* (list)
        """
        distributions,means,stds,moments = self.getSpeciesDistributions(distributions,species)
        self.species_baby_at_division_distributions = distributions
        self.species_baby_at_division_means = means
        self.species_baby_at_division_standard_deviations = stds
        self.species_baby_at_division_moments = moments
       
        
    def setVolumeExtantDistribution(self,distribution):
        """
        Set the volume array
        
        Input:
         - *distribution* (list)
        """
        distribution,mean,std,moments = self.getVolumeDistribution(distribution)
        self.volume_extant_distribution = distribution
        self.volume_extant_mean = mean
        self.volume_extant_standard_deviation = std
        self.volume_extant_moments = moments
        self.HAS_VOLUME_EXTANT = True
        
        
    def setVolumeExtantAtBirthDistribution(self,distribution):
        """
        Set the volume array
        
        Input:
         - *distribution* (list)
        """
        distribution,mean,std,moments = self.getVolumeDistribution(distribution)
        self.volume_extant_at_birth_distribution = distribution
        self.volume_extant_at_birth_mean = mean
        self.volume_extant_at_birth_standard_deviation = std
        self.volume_extant_at_birth_moments = moments              


    def setVolumeExtantAtDivisionDistribution(self,distribution):
        """
        Set the volume array
        
        Input:
         - *distribution* (list)
        """
        distribution,mean,std,moments = self.getVolumeDistribution(distribution)
        self.volume_extant_at_division_distribution = distribution
        self.volume_extant_at_division_mean = mean
        self.volume_extant_at_division_standard_deviation = std
        self.volume_extant_at_division_moments = moments
        

    def setVolumeBabyAtDivisionDistribution(self,distribution):
        """
        Set the volume array
        
        Input:
         - *distribution* (list)
        """
        distribution,mean,std,moments = self.getVolumeDistribution(distribution)
        self.volume_baby_at_division_distribution = distribution
        self.volume_baby_at_division_mean = mean
        self.volume_baby_at_division_standard_deviation = std
        self.volume_extant_at_division_moments = moments        

       
    def setSpeciesAtDivision(self,species,lbls = None): 
        """
        Set the species array
        
        Input:
         - *species* an array of species vs time data
         - *lbls* [default=None] a list of species labels
        """
        self.species_at_division = np.array(species).astype(np.uint32)
        self.HAS_SPECIES_AT_DIVISION = True
        if lbls != None:
            self.species_labels = lbls        
            
    
    def setSpeciesAtBirth(self,species,lbls = None):
        """
        Set the species array
        
        Input:
         - *species* an array of species vs time data
         - *lbls* [default=None] a list of species labels
        """
        self.species_at_birth = np.array(species).astype(np.uint32)
        self.species_at_birth_not_tracked = self.species_at_division - species
        self.HAS_SPECIES_AT_BIRTH = True
        if lbls != None:
            self.species_labels = lbls
            
            
    def setVolumeAtDivision(self,volume):
        """
        set the volume array at division
        
        Input:
         - *volume*
        """
        self.volume_at_division = np.array(volume)
        self.HAS_VOLUME_AT_DIVISION = True
        

    def setVolumeAtBirth(self, volume_selected, volume_notselected):
        """
        set the volume array at birth
        
        Input:
         - *volume*
        """
        self.volume_at_birth = np.array(volume_selected)
        self.volume_at_birth_notselected = np.array(volume_notselected)        
        self.volume_at_birth_all = np.concatenate((self.volume_at_birth, self.volume_at_birth_notselected))
        self.HAS_VOLUME_AT_BIRTH = True    
        
            
    def setXData(self, xdata, lbls=None):
        """
        Sets an array of extra simulation data

        Input:
        - *xdata* an array of xdata vs time
        - *lbls* [default=None] a list of xdata labels
        """
        self.xdata = xdata
        self.HAS_XDATA = True
        if lbls != None:
            self.xdata_labels = lbls
            

    def setSpeciesAtDivisionDistributions(self,distributions,means,stds,moments):
        """
        set species at division distributions for the lineage
        
        Input:
         - *distributions* (list)
         - *means* (dict)
         - *stds* (dict)
         - *moments* (dict) 
        """
        self.species_at_division_distributions = distributions
        self.species_at_division_means = means
        self.species_at_division_standard_deviations = stds
        self.species_at_division_moments = moments
        

    def setSpeciesAtBirthDistributions(self,distributions,means,stds,moments):
        """
        set species at birth distributions for the lineage
        
        Input:
         - *distributions* (list)
         - *means* (dict)
         - *stds* (dict)
         - *moments* (dict)
        """
        self.species_at_birth_distributions = distributions
        self.species_at_birth_means = means
        self.species_at_birth_standard_deviations = stds
        self.species_at_birth_moments = moments
        
    def setVolumeAtDivisionDistribution(self,distributions,means,stds,moments):
        """
        setSpeciesDist stuff for the determination of distributions
        
        Input:
         - *distributions* (list)
         - *means* (dict)
         - *stds* (dict)
         - *moments* (dict)
        """
        self.volume_at_division_distributions = distributions
        self.volume_at_division_means = means
        self.volume_at_division_standard_deviations = stds
        self.volume_at_division_moments = moments
        

    def setVolumeAtBirthDistribution(self,distributions,means,stds,moments):
        """
        setSpeciesDist stuff for the determination of distributions
        
        Input:
         - *distributions* (list)
         - *means* (dict)
         - *stds* (dict)
         - *moments* (dict)
        """
        self.volume_at_birth_distributions = distributions
        self.volume_at_birth_means = means
        self.volume_at_birth_standard_deviations = stds
        self.volume_at_birth_moments = moments
        

    def getTime(self, lbls=False):
        """
        Return the time vector

        Input:
         - *lbls* [default=False] return only the time array or optionally both the time array and time label
        """
        output = None
        if self.HAS_TIME:
            output = self.time.reshape(len(self.time),)
        if not lbls:
            return output
        else:
            return output, [self.time_label]
                  

    def getSpeciesAtDivision(self, lbls=False):
        """
        Return an array of time+species

        Input:
        - *lbls* [default=False] return only the time+species array or optionally both the data array and a list of column label
        """
        output = None
        if self.HAS_SPECIES_AT_DIVISION:
            length = len(self.species_at_division)  # 20-02 prevents unequal sized arrays in the hstack
            output = np.hstack((self.division_times[:length].reshape(len(self.division_times[:length]),1), self.species_at_division))
            labels = [self.time_label]+self.species_labels
        else:
            output = self.division_times
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels
            
            
    def getSpeciesAtBirth(self, lbls=False):
        """
        Return an array of time+species

        Input:
        - *lbls* [default=False] return only the time+species array or optionally both the data array and a list of column label
        """
        output = None
        if self.HAS_SPECIES_AT_BIRTH:
            length = len(self.species_at_birth)
            output = np.hstack((self.division_times[:length].reshape(len(self.division_times[:length]), 1), self.species_at_birth))
            labels = [self.time_label]+self.species_labels
        else:
            output = self.division_times
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels            

    def setDeterministicData(self,time,volume,ages):
        """
        Set the time and volume arrays
        
        Input:
         - *time* (list)
         - *volume* (list)
         - *ages* (list
        """
        self.time_deterministic = np.array(time)
        self.volume_deterministic = np.array(volume)
        self.ages_deterministic = np.array(ages)
        self.HAS_DETERMINISTIC = True     
            

    def getVolumeDeterministic(self): # TODO: move?
        """ Return an array of time+volume """        
        if self.HAS_DETERMINISTIC:            
            output = np.column_stack((self.time_deterministic, self.volume_deterministic))
        else:
            output = self.time_determistic
        return output   
        
            
    def getVolumeAtDivision(self):
        """ Return an array of time+volume """
        output = None
        if self.HAS_VOLUME_AT_DIVISION:            
            output = np.column_stack((self.division_times, self.volume_at_division))          
        else:
            output = self.division_times      
        return output
        
        
    def getVolumeAtBirth(self):
        """ Return an array of time+volume """
        output = None
        if self.HAS_VOLUME_AT_BIRTH:  
            length = len(self.volume_at_birth)   
            output = np.column_stack((self.division_times[:length], self.volume_at_birth))        
        else:
            output = self.division_times
        return output        
        

    def getDataAtTime(self, time):
        """
        Return all data generated at "time"

        Input:
         - *time* the required exact time point
        """        
        t = None
        sp = None        
        temp_t = self.time.reshape(len(self.time),)
        for tt in range(len(temp_t)):
            if temp_t[tt] == time:
                t = tt
                if self.HAS_AGES:
                    sp = self.ages.take([tt], axis=0)                
                break

        output = None
        if t is not None:
            output = np.array([[temp_t[t]]])
            if sp is not None:
                output = np.hstack((output,sp))            
        return output
           

    def getDataInTimeInterval(self, time, bounds=None):
        """
        Returns an array of all data in interval: time-bounds <= time <= time+bounds
        where bound defaults to step size

        Input:
         - *time* the interval midpoint
         - *bounds* [default=None] interval half span defaults to step size
        """
        temp_t = self.time.reshape(len(self.time),)
        if bounds == None:
            bounds = temp_t[1] - temp_t[0]
        c1 = (temp_t >= time-bounds)
        c2 = (temp_t <= time+bounds)
        print('Searching ({0}:{1}:{2})'.format(time-bounds, time, time+bounds))

        t = []
        sp = None
        ra = None
        for tt in range(len(c1)):
            if c1[tt] and c2[tt]:
                t.append(tt)
        output = None
        if len(t) > 0:
            output = self.time.take(t)
            output = output.reshape(len(output),1)
            if self.HAS_AGES and self.HAS_TIME:
                output = np.hstack((output, self.ages.take(t, axis=0)))
        return output
        
        
