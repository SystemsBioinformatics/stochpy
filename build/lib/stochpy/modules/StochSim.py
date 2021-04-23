 #! /usr/bin/env python
"""
Stochastic Simulation Module
============================

;;
Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: Augustus 10, 2015
"""

############################ IMPORTS ################################
from __future__ import division, print_function, absolute_import

import sys,copy,time,os,subprocess,bisect,math,shutil
from io import BytesIO   

try:
    import pickle
except ImportError:
    import cPickle as pickle    

#from . import Analysis
from .StochPyPlot import *
from .StochPyPrint import PrintingFunctions
from .PyscesMiniModel import IntegrationStochasticDataObj,RegularGridDataObj
from ..tools.Progress_bar import Progress_bar
from ..tools.ParseDistributions import ParseDistributions,convertInput2Indices

try: 
    from . import InterfaceCain    
    IS_STOCHPY_CAIN = True
except ImportError: 
    IS_STOCHPY_CAIN = False
    
try: 
    from . import InterfaceStochKit
    from . import PSC2StochkitXML
    InterfaceStochKit.DeleteExistingData()
    IS_STOCHPY_KIT = True
except ImportError:
    IS_STOCHPY_KIT = False
    

class SSASettings():
    """   
    Input:
     - *x_matrix* (array)
     - *timesteps* (integer)
     - *starttime* (float)
     - *endtime* (float)
     - *track_propensities* (boolean)
     - *species_selection* (list)
     - *rate_selection* (list)
     - *last_timepoint* (boolean)
     - *seed* (boolean)
     - *quiet* (boolean)     
    """
    def __init__(self,x_matrix,timesteps,starttime,endtime,track_propensities,species_selection,rate_selection,last_timepoint,seed,quiet):
        self.X_matrix = x_matrix
        self.timesteps = timesteps
        self.starttime = starttime
        self.endtime = endtime 
        self.IsTrackPropensities = track_propensities
        self.species_selection = species_selection        
        self.rate_selection = rate_selection
        self.IsOnlyLastTimepoint = last_timepoint
        self.IsSeed = seed
        self.quiet = quiet

############################ END IMPORTS ############################

class SSA(PlottingFunctions,PrintingFunctions):
    """    
    Input options:
     - *method* [default = 'Direct'], Available methods: 'Direct', 'FirstReactionMethod','TauLeaping','Next Reaction Method'
     - *model_file* [default = 'ImmigrationDeath.psc']
     - *dir* [default = '/home/user/stochpy/pscmodels/ImmigrationDeath.psc']
     - *mode* [default = 'steps'] (string) as simulation for a total number of 'steps' or until a certain end 'time' 
     - *end* [default = 1000] (float) as end of the simulation (number of steps or end time)  
     - *trajectories* [default = 1] (integer)   
     - *IsTrackPropensities* [default = False] (boolean)
     - *IsInteractive* [default = False] (boolean)
     - *IsSeed* [default = False] (boolean)
     - *IsQuiet* [default = True] (boolean)
    """
    def __init__(self,method='direct',model_file='ImmigrationDeath.psc',dir=None,mode='steps',end=1000,trajectories=1,IsTrackPropensities=False,IsInteractive=True,IsSeed=False,IsQuiet=True):
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
        
        self.model_file = model_file
        self.model_dir = dir
        self.output_dir = output_dir
        self.temp_dir = temp_dir    
        self.sim_end = end
        self.sim_mode = mode
        self.sim_trajectories = trajectories    
        self._IsTrackPropensities = IsTrackPropensities
        self._IsSeed = IsSeed
        self._IsQuiet = IsQuiet
        self._IsSimulationDone = False
        self._IsFixedIntervalMethod = False        
        self._IsOnlyLastTimepoint = False
        self._IsTrajectoriesSetByUser = False
        self._MethodSetBy = "DoStochSim"  #08-01-2014
        self._IsModeSetByUser = False
        self._IsEndSetByUser = False
        self._IsDeletePreviousSimulationData = True      
        self._IsCellDivision = False  
        self.HAS_AVERAGE = False
        self.HAS_DELAY_PARAMETERS = False #10-11-2013           
        self.HAS_PUTATIVE_REACTION_TIMES = False 
        self.species_amount_modifications = {}
        self.parameter_value_modifications = {}
        self.Method(method)
        if IsInteractive:          
            try: 
                Analysis.plt.ion()   # Set on interactive pylab environment              
            except Exception as er:
                print(er)
        else:
            try: 
                Analysis.plt.ioff()  # Set on interactive pylab environment
            except Exception as er:
                print(er) 


    def SetSeeding(self,seed=True):
        """
        Input:
         - *seed* [default=True] (boolean)        
        """
        self._IsSeed = seed
        

    def SetQuiet(self,quiet=True):
        """
        Input:
         - *quiet* [default=True] (boolean)
        """
        self._IsQuiet = quiet
                            

    def Method(self,method):
        """       
        Input:
         - *method* (string)    ['direct','frm','nrm','tauleap','DelayedDirect','DelayedNRM','SMM','fSMM']   
         
        Both the SMM and fSMM are single molecule methods and under development
        """
        self._IsTauleaping = False
        self._IsDelayedMethod = False
        self._IsSingleMoleculeMethod = False
        self._IsNRM = False        
                
        method = method.lower()         # Capital insensitive
        if method == 'direct':
            from ..implementations import DirectMethod
            self.sim_method = DirectMethod.DirectMethod
            if not self._IsQuiet:
                print("Info: Direct method is selected to perform stochastic simulations.")
            self.sim_method_name = "Direct"
        elif method in ['firstreactionmethod','frm']:
            from ..implementations import FirstReactionMethod as FRM
            self.sim_method = FRM.FirstReactionMethod
            self.sim_method_name = "FirstReactionMethod"
            if not self._IsQuiet:
                print("Info: First Reaction method is selected to perform stochastic simulations.")
        elif method in ['tauleaping','tau-leaping','tau-leap','tauleap']:
            from ..implementations import TauLeaping
            self.sim_method = TauLeaping.OTL
            self.sim_method_name = "TauLeaping"
            if not self._IsQuiet:
                print("Info: Optimized Tau-Leaping method is selected to perform stochastic simulations.")
                print("Info: User can change the 'epsilon' parameter with DoStochSim(epsilon = 0.01)")
            self._IsTauleaping = True
        elif method in ['nextreactionmethod','nrm']:
            from ..implementations import NextReactionMethod as NRM
            self.sim_method = NRM.NextReactionMethod 
            if not self._IsQuiet:
                print("Info: Next Reaction method is selected to perform stochastic simulations")
            self._IsNRM = True
            self.sim_method_name = "NextReactionMethod"
        elif method in ['delayeddirect','delayed']:
            from ..implementations import DelayedDirectMethod
            self.sim_method = DelayedDirectMethod.DelayedDirectMethod
            if not self._IsQuiet:
                print("Info: Delayed Direct Method is selected to perform delayed stochastic simulations.")
            self.sim_method_name = "DelayedDirectMethod"
            self._IsDelayedMethod = True
        elif method in ['delayednextreactionmethod','delayednrm']:
            from ..implementations import DelayedNRM
            self.sim_method = DelayedNRM.DelayedNRM 
            if not self._IsQuiet:
                print("Info: Delayed Next Reaction Method is selected to perform delayed stochastic simulations.")
            self.sim_method_name = "DelayedNextReactionMethod"
            self._IsDelayedMethod = True
            self._IsNRM = True       
        elif method in ['singlemoleculemethod','smm']:
            from ..implementations import SingleMoleculeMethod
            self.sim_method = SingleMoleculeMethod.SingleMoleculeMethod
            if not self._IsQuiet:
                print("Info: full Single Molecule Method with support for second-order reactions, is selected to perform delayed stochastic simulations.")
                print("*** Warning: The Single Molecule Method assumes that each reaction is described with mass-action kinetics. Use the fSMM if there are non-linear rate equations ***")           
            self.sim_method_name = "SingleMoleculeMethod"            
            self._IsSingleMoleculeMethod = True            
        elif method in ['fastsinglemoleculemethod','fsmm']: 
            from ..implementations import FastSingleMoleculeMethod
            self.sim_method = FastSingleMoleculeMethod.FastSingleMoleculeMethod 
            if not self._IsQuiet:
                print("Info: fast Single Molecule Method is selected to perform delayed stochastic simulations.")            
            self.sim_method_name = "FastSingleMoleculeMethod"
            self._IsNRM = True
            self._IsSingleMoleculeMethod = True
        else:
            print("*** WARNING ***: Only valid options are: 'Direct', 'FirstReactionMethod', 'NextReactionMethod','TauLeap' or 'DelayedDirect', 'DelayedNRM', 'SMM', 'fSMM'.")
            print("Info: By default, the Direct method is selected")
            from ..implementations import DirectMethod
            self.sim_method = DirectMethod.DirectMethod
            self.sim_method_name = "Direct"   
        
        self.SSA = self.sim_method(self.model_file,self.model_dir,self._IsQuiet)
        for s in self.species_amount_modifications:
            self.ChangeInitialSpeciesCopyNumber(s,self.species_amount_modifications[s])            
        for p in self.parameter_value_modifications: 
            self.ChangeParameter(p,self.parameter_value_modifications[p])
        
        self.data_stochsim = IntegrationStochasticDataObj()   
        self.data_stochsim_grid = RegularGridDataObj()   
        self._IsSimulationDone = False    
        self.HAS_AVERAGE = False   
        self._MethodSetBy = "User"        
      

    def Timesteps(self,s):
        """       
        Set the number of time steps to be generated for each trajectory
        
        Input:
         - *s* (integer)
        """              
        try:
            self.sim_end  = abs(int(s))
            self.sim_mode = 'steps'  
            if not self._IsQuiet:                
                print("Info: The number of time steps is: {0:d}".format(self.sim_end))
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
            self.sim_end = abs(float(t))
            self.sim_mode = 'time'            
            if not self._IsQuiet:           
                print("Info: The simulation end time is: {0}".format(self.sim_end))
            self._IsEndSetByUser = True
            self._IsModeSetByUser = True
        except ValueError:
            raise ValueError("The end time must be an integer or a float")


    def Trajectories(self,n):
        """
        Set the number of trajectories to be generated
        
        Input:
         - *n* (integer)
        """
        if isinstance(n,int):  
            self.sim_trajectories = n
            self._IsTrajectoriesSetByUser = True
        else:
            raise TypeError("Integer argument expected, got {0}".format(type(n)))     


    def Reload(self):
        """ Reload the entire model again. Useful if the model file has changed """
        self.SSA.Parse(self.model_file,self.model_dir,IsTauleaping = self._IsTauleaping, IsNRM = self._IsNRM,IsDelayed = self._IsDelayedMethod, IsSMM = self._IsSingleMoleculeMethod,IsQuiet=self._IsQuiet)              
        self._IsSimulationDone = False   
        self.HAS_AVERAGE = False 
        self.data_stochsim = IntegrationStochasticDataObj()       
        self.data_stochsim_grid = RegularGridDataObj()  
        
        self.species_amount_modifications = {}
        self.parameter_value_modifications = {}
        self.sim_trajectories_done = None
        self._current_trajectory = None
        self.sim_species_tracked = None
        self.settings = None


    def Model(self,model_file,dir=None):
        """        
        Select model for simulation
              
        Input:
         - *model_file* (string)
         - *dir* [default = None] (string)
        """      
        if self.HAS_DELAY_PARAMETERS:
            self.HAS_DELAY_PARAMETERS = False
            if not self._IsQuiet:      
                print('Info: Please mind that the delay parameters have to be set again.')
            self.delay_distr_parameters = None
            self.delay_distributions = None
            self.delayed_nonconsuming = None
            self.delayed_consuming = None
        if self.HAS_PUTATIVE_REACTION_TIMES:
            self.HAS_PUTATIVE_REACTION_TIMES = False
            if not self._IsQuiet:
                print('Info: Please mind that the waiting time distributions have to be set again.')
            self._putative_reaction_times_distr_parameters = None
            self._putative_reaction_times = None
        self.model_file = model_file
        if dir != None:
            self.model_dir = dir
        self.Reload()


    def Mode(self,sim_mode='steps'):
        """      
        Run a stochastic simulation for until `end` is reached. This can be either time steps or end time (which could be a *HUGE* number of steps).

        Input:
         - *sim_mode* [default = 'steps'] (string) 'time' or 'steps'      
        """
        self.sim_mode = sim_mode.lower()
        if self.sim_mode.lower() not in ['steps','time']:
            print("*** WARNING ***: Mode '{0}' not recognized using: 'steps'".format(sim_mode))
            self.sim_mode = 'steps'        
        self._IsModeSetByUser = True


    def GetTrajectoryData(self,n=1):
        """        
        Switch to another trajectory, by default, the last trajectory is directly accessible      
        
        Input:
         - *n* [default = 1] (integer)
        """ 
        if isinstance(n,int):
            try:      
                file_in = open(os.path.join(self.temp_dir,'{0}{1:d}.dat'.format(self.model_file,n)),'rb') # rb necessary for Python 3.x
                self.data_stochsim = pickle.load(file_in)
                file_in.close()
            except IOError:
                raise IOError("Trajectory {0:d} does not exist".format(n))                
        else:
            raise TypeError("Integer argument expected, got float")


    def DumpTrajectoryData(self,n):
        """
        Input:
         - *n* (integer)
        """ 
        if isinstance(n,int):    
            if n == self.data_stochsim.simulation_trajectory:                        
                filename_out = os.path.join(self.temp_dir,'{0}{1:d}.dat'.format(self.model_file,n))                          
                f = open(filename_out,'wb')
                pickle.dump(self.data_stochsim,f)         
                f.close()
            else:
                print("Error: Trajectory {0} is currently not selected/ does not exist".format(n))
        else:
            raise TypeError("Integer argument expected, got float")
            

    def ChangeParameter(self,parameter,value):
        """        
        Change parameter value   
        
        Input:
         - *parameter* (string)
         - *value* (float)
        """
        IsKeyError = False
        if isinstance(parameter,str) and type(value) in [int,float,np.float64,np.float32]:   
            try:
                self.SSA.parse.Mod.__pDict__[parameter]['initial'] = float(value)               
            except KeyError:              
                print("Parameters are: {0}".format(list(self.SSA.parse.Mod.__pDict__)))
                IsKeyError = True
            if not IsKeyError: 
                self.SSA.parse.BuildReactions()
                self.SSA.propensities = copy.deepcopy(self.SSA.parse.propensities)           
        else:
            print("Error: arguments parameter = string and value = float")
        
        self.parameter_value_modifications[parameter] = value


    def ChangeInitialSpeciesAmount(self,species,value):
        """
        Change initial species copy number
        
        Input:
         - *species* (string)
         - *value* (float)
        """
        print('\n*** Deprecation warning ***:\nChangeInitialSpeciesAmount() is being replaced with ChangeInitialSpeciesCopyNumber()')
        self.ChangeInitialSpeciesCopyNumber(species,value)        
 

    def ChangeInitialSpeciesCopyNumber(self,species,value):
        """
        Change initial species copy number
        
        Input:
         - *species* (string)
         - *value* (float)
        """
        IsKeyError = False
        if isinstance(species,str) and type(value) in [int,float,np.float64,np.float32]:   
            try:
                self.SSA.parse.Mod.__sDict__[species]['initial'] = float(value)
            except KeyError:     
                print("Species are: {0}".format(list(self.SSA.parse.Mod.__sDict__)))
                IsKeyError = True              
            if not IsKeyError:
                if self.SSA.parse.Mod.__sDict__[species]['fixed']: # rebuild reactions and propensities
                    self.SSA.parse.BuildReactions()
                    self.SSA.propensities = copy.deepcopy(self.SSA.parse.propensities)
                self.SSA.parse.BuildX()     
                self.SSA.X_matrixinit = copy.deepcopy(self.SSA.parse.X_matrix.transpose()[0])
        else:
            print("Error: species argument must be a string and value argument must be a float or integer") 
        
        self.species_amount_modifications[species] = value          
                

    def DoStochKitStochSim(self,endtime=100,frames=10000,trajectories=False,IsTrackPropensities=False,customized_reactions=None,solver=None,keep_stats = False,keep_histograms = False,quiet = False):
        """
        Do Stochastic simulations with StochKit in StochPy
        Make sure that the input file contains no net stoichiometric coefficients
        
        Input:
         - *endtime* [default = 100] (float)
         - *frames* [default = 10000] (integer)
         - *trajectories* [default = False] (integer)
         - *IsTrackPropensities* [default = False] (boolean)
         - *customized_reactions* [default=None] (list of strings)
         - *solver* [default = None] (string)
         - *keep_states* [default = False] (boolean)
         - *keep_histograms* [default = False) (boolean)
         - *quiet* [default = False] (boolean) suppress StochPy print statements
        """
        if self._IsQuiet:
            quiet = True        
        
        assert IS_STOCHPY_KIT,"StochKit and/or InterfaceStochKit is not installed or the directories in InterfaceStochKit.ini are incorrect"        
        if IS_STOCHPY_KIT:               
            self.DeleteTempfiles() # Delete '.dat' files
            if not quiet: print("*** WARNING ***: Do not use net stoichiometric coefficients for fixed-interval output solvers.")
            if isinstance(frames,int):            
                pass
            elif type(frames) in [float,np.float64,np.float32]:
                if not quiet: print("*** WARNING ***: 'frames' must be an integer rather than a float; float {0} is rounded to {1:d}".format(frames,int(frames)))
                frames = int(frames) 
            else:
                print("Error: 'frames' must be an integer")
                sys.exit()   
            
            self.data_stochsim = IntegrationStochasticDataObj()
            self.data_stochsim_grid = RegularGridDataObj()
            
            self._IsOnlyLastTimepoint = False
            self._IsFixedIntervalMethod = True
            self._IsSimulationDone = False
            self._IsTrackPropensities = IsTrackPropensities
            if self._IsTrackPropensities:
                self.sim_rates_tracked = copy.copy(self.SSA.rate_names) 
            self.HAS_AVERAGE = False          
            if trajectories != False:
                self.Trajectories(trajectories)             
                self._IsTrajectoriesSetByUser = False
            elif trajectories == False and self.sim_trajectories != 1 and not self._IsTrajectoriesSetByUser:
                    self.Trajectories(1)
            if customized_reactions != None:
                for r_id in customized_reactions:                   
                    self.SSA.parse.Mod.__nDict__[r_id]['stochtype'] = 'customized'                    
                    
            assert not self.SSA.parse.Mod.__HAS_ASSIGNMENTS__, "StochKit solvers do not support assignments. Use the StochPy solvers DoStochSim()"          
   
            if solver == None:                  
                solver = InterfaceStochKit.STOCHKIT_SOLVER # use the default solver
            
            t1 = time.time()
            stochkit_model_filename = self.model_file.split('.')[0]+'_stochkit.xml'
            doc = PSC2StochkitXML.xml_createStochKitDoc(self.SSA)
            PSC2StochkitXML.xml_viewStochKitXML(doc,fname=os.path.join(self.temp_dir,stochkit_model_filename))          
            stochkit_model_filename_path = os.path.join(self.temp_dir, stochkit_model_filename)            
            stochkit_keep_stats = keep_stats
            stochkit_keep_histograms = keep_histograms
            stochkit_keep_trajectories = True
            stochkit_cmd = "-m {0} -t {1} -r {2} -i {3} --label -f".format(stochkit_model_filename_path,endtime,self.sim_trajectories ,frames)
            if not stochkit_keep_stats:
                stochkit_cmd += " --no-stats"
            if stochkit_keep_histograms:
                stochkit_cmd += " --keep-histograms"
            if stochkit_keep_trajectories:
                stochkit_cmd += " --keep-trajectories"
            stochkit_cmd += " --out-dir {0}".format(os.path.join(InterfaceStochKit.STOCHKIT_WORK_DIR,stochkit_model_filename))
            #print stochkit_cmd
            if not quiet:            
                if self.sim_trajectories == 1:                  
                    print("Info: {0:d} trajectory is generated until t = {1} with {2:d} frames".format(self.sim_trajectories,endtime,frames))
                else:                  
                    print("Info: {0:d} trajectories are generated until t = {1} with {2:d} frames".format(self.sim_trajectories,endtime,frames))    
            try: 
                solver_path = os.path.join(InterfaceStochKit.STOCHKIT_SOLVER_DIR, solver)                
                #rcode = subprocess.call([solver_path, stochkit_cmd]) 
                _string ='{0:s} {1:s}'.format(solver_path,stochkit_cmd)
                rcode = subprocess.call(_string.split())
                IsSimulation = True
            except Exception as er:
                print(er)
                print(solver_path)
                IsSimulation = False
                                     
            if IsSimulation:    
                t2 = time.time()    
                self.simulation_time = t2-t1    
                if not quiet:
                    print("Info: Simulation time including compiling {0:1.5f}".format(self.simulation_time))
                    if self._IsTrackPropensities:
                        print("Info: Parsing data to StochPy and calculating propensities and distributions ...")
                    else:
                        print("Info: Parsing data to StochPy and calculating distributions ...")
                self.data_stochsim = InterfaceStochKit.GetStochKitOutput(stochkit_model_filename,self.SSA,self.model_file,endtime,self.sim_trajectories,frames,self._IsTrackPropensities)              
                self.sim_species_tracked = copy.copy(self.SSA.species_names)    
                self.sim_trajectories_done = copy.copy(self.sim_trajectories)
                try:
                    self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.SSA.rate_names,self.plot.plotnum,quiet)
                except:
                    self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.SSA.rate_names,quiet = quiet)
                self._IsSimulationDone = True
                if not quiet:
                   print("Info: Data successfully parsed into StochPy")
       

    def DoCainStochSim(self,endtime=100,frames=10000,trajectories=False,solver="HomogeneousDirect2DSearch",IsTrackPropensities=False,quiet = False):
        """      
        Use Cain implementations for fixed-interval output stochastic simulations (www.cacr.caltech.edu/~sean/cain/DeveloperFile)
        Make sure that the input file contains no net stoichiometric coefficients
        
        Input:
         - *endtime* [default = 100](float)
         - *frames* [default = 10000] (integer)
         - *trajectories* [default = False] (integer)
         - *solver* [default = 'HomogeneousDirect2DSearch'] (string)
         - *IsTrackPropensities* [default = False] (boolean)
         - *quiet* [default = False] (boolean) suppress StochPy print statements
        """   
        if self._IsQuiet:
            quiet = True     
        assert IS_STOCHPY_CAIN, "InterfaceCain is not installed"      
        if IS_STOCHPY_CAIN:
            if not quiet: print("*** WARNING ***: Only mass-action kinetics can be correctly parsed by the Cain solvers")          
            self.DeleteTempfiles() # Delete '.dat' files
            if not quiet: print("*** WARNING ***: Do not use net stoichiometric coefficients for fixed-interval output solvers.")
            if isinstance(frames,int):
                pass
            elif type(frames) in [float,np.float64,np.float32]:
                if not quiet: print("*** WARNING ***: 'frames' must be an integer rather than a float; float {0} is rounded to {1:d}".format(frames,int(frames)))
                frames = int(frames) 
            else:
                print("Error: 'frames' must be an integer")
                sys.exit()             

            self.data_stochsim = IntegrationStochasticDataObj()
            self.data_stochsim_grid = RegularGridDataObj()         
        
            self._IsOnlyLastTimepoint = False
            self._IsFixedIntervalMethod = True
            self._IsTrackPropensities = IsTrackPropensities
            if self._IsTrackPropensities:
                self.sim_rates_tracked = copy.copy(self.SSA.rate_names)    
            self.HAS_AVERAGE = False  
            self._IsSimulationDone = False          
            if trajectories != False:
                self.Trajectories(trajectories)
                self._IsTrajectoriesSetByUser = False
            elif trajectories == False and self.sim_trajectories != 1 and not self._IsTrajectoriesSetByUser:
                self.Trajectories(1)
                      
            assert not self.SSA.parse.Mod.__HAS_EVENTS__, "Cain solvers do not support events. Use the StochPy solvers DoStochSim()"               
            assert not self.SSA.parse.Mod.__HAS_ASSIGNMENTS__, "Cain solvers do not support assignments. Use the StochPy solvers DoStochSim()"
            ### Parse model to CAIN ###
            mersenne_twister_data = InterfaceCain.getCainInputfile(self.SSA,endtime,frames,self.sim_trajectories)            
            cain_cmd_filename = 'cain_in.txt'
            cmd_file = open(os.path.join(self.temp_dir, cain_cmd_filename), 'rb')
            cain_cmd = cmd_file.read()
            cmd_file.close()
            ###########################
            
            try:
                if (os.sys.platform == 'win32') and (not solver.endswith('.exe')):
                    solver = solver.split('.')[0] + '.exe'                 
                solver_path = os.path.join(InterfaceCain.CAIN_SOLVER_PATH,solver)  
                proc = subprocess.Popen(os.path.join(InterfaceCain.CAIN_SOLVER_PATH,solver),stdin=subprocess.PIPE,stdout=subprocess.PIPE) 
                IsFoundSolver = True
            except Exception as er:
                print(er)
                print(solver_path)
                IsFoundSolver = False    
            
            if IsFoundSolver:
                if not quiet:
                    if self.sim_trajectories == 1:                  
                        print("Info: {0:d} trajectory is generated until t = {1} with {2:d} frames".format(self.sim_trajectories,endtime,frames))
                    else:                  
                        print("Info: {0:d} trajectories are generated until t = {1} with {2:d} frames".format(self.sim_trajectories,endtime,frames))
                t1 = time.time()     
                stdout_values = proc.communicate(cain_cmd)[0]              
                t2 = time.time() 
                self.simulation_time = t2-t1            
                if not quiet:
                    print("Info: Simulation time {0:1.5f}".format(self.simulation_time))
                    if self._IsTrackPropensities:
                        print("Info: Parsing data to StochPy and calculating propensities and distributions ...")
                    else:
                        print("Info: Parsing data to StochPy and calculating distributions ...")                    
                
                self.data_stochsim = InterfaceCain.getCainOutput2StochPy(BytesIO(stdout_values).readlines() ,mersenne_twister_data,self.SSA,self.model_file,endtime,self.sim_trajectories,frames,self._IsTrackPropensities)              
                self.sim_species_tracked = copy.copy(self.SSA.species_names)    
                self.sim_trajectories_done = copy.copy(self.sim_trajectories)
                try: 
                    self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.SSA.rate_names,self.plot.plotnum)
                except:
                    self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.SSA.rate_names)
                self._IsSimulationDone = True
                if not quiet: print("Info: Data successfully parsed into StochPy")
    

    def DoStochSim(self,end=False,mode=False,method=False,trajectories=False,epsilon = 0.03,IsTrackPropensities=False, rate_selection = None, species_selection = None,IsOnlyLastTimepoint = False,critical_reactions=[],reaction_orders = False,species_HORs = False,species_max_influence = False,quiet = False):
        """
        Run a stochastic simulation for until `end` is reached. This can be either time steps or end time (which could be a *HUGE* number of steps).

        Input:
         - *end* [default=1000] (float) simulation end (steps or time)
         - *mode* [default='steps'] (string) simulation mode, can be one of: ['steps','time']
         - *method* [default='Direct'] (string) stochastic algorithm ['Direct', 'FRM', 'NRM', 'TauLeap']
         - *trajectories* [default = 1] (integer)
         - *epsilon* [default = 0.03] (float) parameter for the tau-leap method
         - *IsTrackPropensities* [default = False]
         - *rate_selection* [default = None] (list) of names of rates to store. This saves memory space and prevents Memory Errors when propensities propensities are tracked
         - *species_selection* [default = None] (list) of names of species to store. This saves memory space and prevents Memory Errors (occurring at ~15 species).
         - *IsOnlyLastTimepoint* [default = False] (boolean)
         - *critical_reactions* [default = [] ] (list) ONLY for the tau-leaping method where the user can pre-define reactions that are "critical". Critical reactions can fire only once per time step.
         - *reaction_orders* [default = [] (list) ONLY for the tau-leaping method 
         - *species_HORs* [default = []  (list) ONLY for the tau-leaping method 
         - *species_max_influence* [default = []]  (list) ONLY for the tau-leaping method 
         - *quiet* [default = False] suppress print statements
        """
        if self._IsQuiet:
            quiet = True
        if species_selection and isinstance(species_selection,str):   
            species_selection = [species_selection]
        if species_selection and isinstance(species_selection,list): 
            for s_id in species_selection:
                assert s_id in self.SSA.species_names, "Species {0} is not in the model or species selection".format(s_id)
        
        self._IsTrackPropensities = IsTrackPropensities
        if rate_selection and isinstance(rate_selection,str):   
            rate_selection = [rate_selection]
            self._IsTrackPropensities = True
        if rate_selection and isinstance(rate_selection,list): 
            for r_id in rate_selection:
                assert r_id in self.SSA.rate_names, "Reaction {0} is not in the model or reaction selection".format(r_id)
            self._IsTrackPropensities = True
               
        self._IsOnlyLastTimepoint = IsOnlyLastTimepoint
        if mode != False:
            self.Mode(sim_mode = mode) 
            self._IsModeSetByUser = False
        elif mode == False and self.sim_mode != 'steps' and not self._IsModeSetByUser:
            self.Mode('steps')
        if end != False:    
            if type(end) in [int,float,np.float64,np.float32]: 
                self.sim_end = end   
            else:
                print("*** WARNING ***: 'end' should be an integer or float\n 1000 is used by default")
                self.sim_end = 1000   
            self._IsEndSetByUser=False
        elif end == False and self.sim_end != 1000 and not self._IsEndSetByUser:
            self.sim_end = 1000

        self.data_stochsim = IntegrationStochasticDataObj()
        self.data_stochsim_grid = RegularGridDataObj()  
        
        if method != False: 
            self.Method(method)
            self._MethodSetBy = "DoStochSim"
        elif method == False and self.sim_method_name != "Direct" and self._MethodSetBy == "DoStochSim":
            self.Method("Direct")
        
        # If DoStochSim not called from DoDelayedStochSim, the method should never be delayed.
        if (self._IsDelayedMethod or self._IsSingleMoleculeMethod) and self._MethodSetBy != "Script": # 08-01-2014
            print("*** WARNING ***: ({0:s}) was selected. Switching to the direct method.".format(self.sim_method_name))
            self.Method('Direct')
            self._MethodSetBy = "DoStochSim"
        
        if reaction_orders != False:
            if not quiet:
                print("Info: reaction orders {0} replaced with {1}".format(self.SSA.parse.reaction_orders, reaction_orders))
            self.SSA.parse.reaction_orders = reaction_orders
        if species_HORs != False:
            self.SSA.parse.species_HORs = species_HORs
            if not quiet:                
                print("Info: species HORs {0} replaced with {1}".format(self.SSA.parse.species_HORs, species_HORs))
        if species_max_influence != False:
            self.SSA.parse.species_max_influence = species_max_influence
            if not quiet: 
                print("Info: species max influence {0} replaced with {1}".format(self.SSA.parse.species_max_influence, species_max_influence))
                
        if trajectories != False: 
            self.Trajectories(trajectories)                        
            self._IsTrajectoriesSetByUser = False
        elif trajectories == False and self.sim_trajectories != 1 and not self._IsTrajectoriesSetByUser:
            self.Trajectories(1)        
        
        self._IsFixedIntervalMethod = False      
        self.HAS_AVERAGE = False      
        if self._IsDeletePreviousSimulationData:
           self.DeleteTempfiles()  # Delete '.dat' files
      
        if not quiet: 
            if self.sim_trajectories == 1:
                print("Info: 1 trajectory is generated")
            else:      
                print("Info: {0:d} trajectories are generated".format(self.sim_trajectories))
                print("Info: Time simulation output of the trajectories is stored at {0:s} in directory: {1:s}".format(self.model_file[:-4]+'(trajectory).dat',self.temp_dir))
                        
        progressBar = Progress_bar(cycles_total = self.sim_trajectories, done_msg = 'time')  ##Progress bar addition## Shows Simulation time afterwards
        for self._current_trajectory in range(1,self.sim_trajectories+1):
            if self.sim_trajectories > 1:
               IsStatusBar = False
            else:
               IsStatusBar = True
               t1 = time.time() 
            
            if self.sim_mode.lower() == 'time':
                self.settings = SSASettings(x_matrix=self.SSA.X_matrixinit,timesteps=10**50,starttime=0,endtime=self.sim_end, track_propensities=self._IsTrackPropensities, species_selection=species_selection,rate_selection = rate_selection,last_timepoint=IsOnlyLastTimepoint,seed = self._IsSeed,quiet = quiet)              
            elif self.sim_mode.lower() == 'steps':
                self.settings = SSASettings(x_matrix=self.SSA.X_matrixinit,timesteps=self.sim_end,starttime=0,endtime=10**50, track_propensities=self._IsTrackPropensities, species_selection=species_selection,rate_selection = rate_selection,last_timepoint=IsOnlyLastTimepoint,seed = self._IsSeed,quiet = quiet)               
            else:
                print("*** WARNING ***: Simulation mode should be 'time' or 'steps'. Steps is done by default")
                self.settings = SSASettings(x_matrix=self.SSA.X_matrixinit,timesteps=self.sim_end,starttime=0,endtime=10**50, track_propensities=self._IsTrackPropensities, species_selection=species_selection,rate_selection = rate_selection,last_timepoint=IsOnlyLastTimepoint,seed = self._IsSeed,quiet = quiet) 
                
            if self.sim_method_name.lower() == "tauleaping":    
                self.SSA.Execute(self.settings,IsStatusBar,epsilon,critical_reactions)
            else:
                self.SSA.Execute(self.settings,IsStatusBar)

            self.data_stochsim = IntegrationStochasticDataObj()         
                  
            self.FillDataStochsim()          
            if self.sim_trajectories == 1 and not quiet:                 
                print("Info: Number of time steps {0:d} End time {1}".format(self.SSA.timestep,self.SSA.sim_t))
            elif self.sim_trajectories > 1:                
                self.DumpTrajectoryData(self._current_trajectory)
                progressBar.update(quiet=quiet) #Progress bar addition, only for multiple trajectories
               
        self._IsSimulationDone = True      
        self.sim_trajectories_done = copy.copy(self.sim_trajectories)
        try: 
            self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.sim_rates_tracked,self.plot.plotnum,quiet)
        except:
            self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.sim_rates_tracked,quiet=quiet)
        if IsStatusBar:
            t2 = time.time()          
            self.simulation_time = t2-t1
            if not quiet: 
                print("Info: Simulation time {0:1.5f}".format(self.simulation_time))
        else:
            self.simulation_time = progressBar.end_time - progressBar.t1            
        if IsOnlyLastTimepoint and not quiet: 
            print('Info: not enough data points (are stored) to determine statistics.')        
            

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
        # Parse the distribution dictionary to two dictionaries. These have reaction indices as keys and values are resp. distribution functions and parameter lists.
        self.delay_distributions, self.delay_distr_parameters = ParseDistributions(delay_distributions, self.SSA.rate_names)  # use all rate names, because this is done before the simulation starts
        delayed_reactions = list(self.delay_distributions)        
        if nonconsuming_reactions == None:              # ==None to recognize 0 as reaction index.
            self.delayed_nonconsuming = []
            self.delayed_consuming = delayed_reactions  # All reactions consuming            
        else:                                           # Nonconsuming reactions supplied
            self.delayed_nonconsuming = convertInput2Indices(nonconsuming_reactions, self.SSA.rate_names)           
            self.delayed_nonconsuming = [r for r in self.delayed_nonconsuming if r in delayed_reactions] # Selects nonconsuming reactions that have a delay distribution.
            self.delayed_consuming = list(set(delayed_reactions) - set(self.delayed_nonconsuming)) # Selects consuming reactions by: all_reaction - nonconsuming     
                                        
        self.HAS_DELAY_PARAMETERS = True


    def DoDelayedStochSim(self, end=False, mode=False, method=False, trajectories=False, IsTrackPropensities=False, rate_selection = None, species_selection = None, IsOnlyLastTimepoint = False,quiet=False):
        """
        Run a stochastic simulation with delayed reactions until `end` is reached. This can be either time steps or end time (which could be a *HUGE* number of steps).

        Input:
         - *end* [default=1000] simulation end (steps or time)
         - *mode* [default='steps'] simulation mode, can be one of: ['steps','time']
         - *method* [default='DelayedDirect'] stochastic algorithm ['DelayedDirect', 'DelayedNRM']
         - *trajectories* [default = 1]
         - *IsTrackPropensities* [default = False]
         - *rate_selection* [default = None] List of names of rates to store. This saves memory space and prevents Memory Errors when propensities propensities are tracked
         - *species_selection* [default = None] List of names of species to store. This saves memory space and prevents Memory Errors (occurring at ~15 species).
         - *IsOnlyLastTimepoint* [default = False]
         - *quiet* [default = False] suppress print statements
        """
        if self._IsQuiet:
            quiet = True   
        if method != False: 
            self.Method(method)
            self._MethodSetBy = "DoStochSim"
        elif self._MethodSetBy == "DoStochSim" and self.sim_method_name != "DelayedDirectMethod":
            self.Method("DelayedDirect")  # Default
        
        if not self._IsDelayedMethod:
            print("*** WARNING ***: an invalid method ({0}) was selected. Switching to the Delayed Direct Method.".format(self.sim_method_name))
            self.Method('DelayedDirect')  # = Default delayed method      
        
        if self.HAS_DELAY_PARAMETERS:
            # Pass delay parameters to delayed SSA implementation.
            self.SSA.distr_functions        = copy.copy(self.delay_distributions)
            self.SSA.distr_parameters       = copy.copy(self.delay_distr_parameters)
            self.SSA.reactions_Consuming    = copy.copy(self.delayed_consuming)
            self.SSA.reactions_NonConsuming = copy.copy(self.delayed_nonconsuming)
        else:
            raise AttributeError("No delay parameters have been set for the model '{0:s}'. First use the function .SetDelayParameters().".format(self.model_file)) #7-1-2014 exit if no delay parameters
        
        # Specify that delayed method is set by this script. Prevents DoStochSim to select other method.
        temp_MethodSetBy = self._MethodSetBy #Either "User" or "DoStochSim"
        self._MethodSetBy = "Script"    
        
        #Perform Stochastic Simulation with given settings and method selected in
        self.DoStochSim(end=end, mode=mode, method=False, trajectories=trajectories, IsTrackPropensities=IsTrackPropensities, rate_selection = rate_selection, species_selection = species_selection, IsOnlyLastTimepoint = IsOnlyLastTimepoint,quiet=quiet)
        
        # Reset to original value
        self._MethodSetBy = temp_MethodSetBy
        if IsOnlyLastTimepoint and not quiet:
            print('Info: not enough data points (are stored) to determine statistics.')    
    

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
        self._putative_reaction_times, self._putative_reaction_times_distr_parameters = ParseDistributions(distributions, self.SSA.rate_names) 
        self.HAS_PUTATIVE_REACTION_TIMES = True
    

    def DoSingleMoleculeStochSim(self, end=False, mode=False, method=False, trajectories=False, species_selection = None, IsOnlyLastTimepoint = False,quiet=False):
        """     
        Run a single molecule stochastic simulation until `end` is reached. This can be either time steps or end time (which could be a *HUGE* number of steps).

        Input (similar to .DoStochSim()):
         - *end* [default=1000] simulation end (steps or time)
         - *mode* [default='steps'] simulation mode, can be one of: ['steps','time']
         - *method* [default='SingleMoleculeMethod'] stochastic algorithm, can be one of: ['SMM',fSMM']
         - *trajectories* [default = 1]    
         - *species_selection* [default = None] (list) of names of species to store. This saves memory space and prevents Memory Errors.
         - *IsOnlyLastTimepoint* [default = False]
         - *quiet* [default = False] suppress print statements  
        """  
        if self._IsQuiet:
            quiet = True      
        # Check whether waiting time distributions were given.
        if not self.HAS_PUTATIVE_REACTION_TIMES:
            raise Warning("No distributions have been set for the model '{0:s}'. First use the function .SetPutativeReactionTimes().".format(self.model_file))
        
        if method != False: 
            self.Method(method)
            self._MethodSetBy = "DoStochSim"
        elif self._MethodSetBy == "DoStochSim" and self.sim_method_name != "FastSingleMoleculeMethod":
            self.Method("fSMM") # Default method
        
        if not self._IsSingleMoleculeMethod:
            print("*** WARNING ***: an invalid method ({0}) was selected. Switching to the fast Single Molecule Method.".format(self.sim_method_name))
            self.Method('fSMM')
        
        if self.sim_method_name == "FastSingleMoleculeMethod" and 2 in [self.SSA.parse.reaction_orders[j] for j in self._putative_reaction_times]:
            print("*** WARNING ***: Second-order reactions are not supported by the fSMM. Switching to the SMM.")
            self.Method('SMM') 

        # Pass delay parameters to SingleMolecule SSA implementation.
        self.SSA.distr_functions  = copy.copy(self._putative_reaction_times)
        self.SSA.distr_parameters = copy.copy(self._putative_reaction_times_distr_parameters)                    

        #If Single Molecule Method, set exponential distributions to reactions not specified in self._putative_reaction_times.
        if self.sim_method_name == 'SingleMoleculeMethod':# and len(self.SSA.distr_functions) < self.SSA.n_reactions:
            self.SSA.auto_exponential_reactions = []
            #print("\nUsing exponential waiting time distributions for:")          
            for j in range(self.SSA.n_reactions):
                if j not in self.SSA.distr_functions:                       # Don't replace already assigned distributions.
                    self.SSA.distr_functions[j] = np.random.exponential
                    self.SSA.distr_parameters[j] = np.nan                   # 31-03-2014 To be specified at start of simulation (self.SSA.EvaluatePropensities)
                    self.SSA.auto_exponential_reactions.append(j)           # 31-03-2014
       
        # Specify that delayed method is set by the script here. Prevents DoStochSim to select other method.
        temp_MethodSetBy = self._MethodSetBy # Either "User" or "DoStochSim"
        self._MethodSetBy = "Script"           
        
        if self._IsTrackPropensities:
            print("*** WARNING ***: Propensities cannot be tracked with the single molecule method")
            self._IsTrackPropensities = False
        
        self.DoStochSim(end=end, mode=mode, method=False, trajectories=trajectories, IsTrackPropensities=False, species_selection = species_selection, IsOnlyLastTimepoint = IsOnlyLastTimepoint,quiet=quiet)
        self._MethodSetBy = "DoStochSim"     # RESET
        
        # Reset to original value
        self._MethodSetBy = temp_MethodSetBy
        if IsOnlyLastTimepoint and not quiet:
            print('Info: not enough data points (are stored) to determine statistics.')


    def DoCompleteStochSim(self, error = 0.001, size=100000,IsTrackPropensities=False, rate_selection=None, species_selection = None,quiet=False):
        """
        Do a stochastic simulation until the first four moments converge (in development, beta-status)
        
        Input:
         - *error* maximal allowed error [default = 0.001]
         - *size* (integer) number of steps before checking the first four moments [default = 100000]
         - *IsTrackPropensities* [default = False]
         - *rate_selection* [default = None] (list) of names of rates to store (saves memory space and prevents Memory Errors when propensities propensities are tracked)
         - *species_selection* [default = None] (list) of names of species to store  (saves memory space and prevents Memory Errors (occurring at ~15 species))
         - *quiet* [default = False] suppress print statements
        """
        if self._IsQuiet:
            quiet = True   
        if species_selection and isinstance(species_selection,str):   
            species_selection = [species_selection]
        if species_selection and isinstance(species_selection,list): 
            for s_id in species_selection:
                assert s_id in self.SSA.species_names, "Species {0} is not in the model or species selection".format(s_id)
        
        self._IsTrackPropensities = IsTrackPropensities       
        if rate_selection and isinstance(rate_selection,str):   
            rate_selection = [rate_selection]
            self._IsTrackPropensities = True
        if rate_selection and isinstance(rate_selection,list): 
            for r_id in rate_selection:
                assert r_id in self.SSA.rate_names, "Reaction {0} is not in the model or reaction selection".format(r_id)
            self._IsTrackPropensities = True

        self.Trajectories(1)                
        self._IsFixedIntervalMethod = False
        self.HAS_AVERAGE = False      
        self.DeleteTempfiles()    # Delete '.dat' files   
        
        self.data_stochsim = IntegrationStochasticDataObj()
        self.data_stochsim_grid = RegularGridDataObj()          
        
        self._current_trajectory = 1
        t1 = time.time()
        self.settings = SSASettings(x_matrix=self.SSA.X_matrixinit,timesteps=size,starttime=0,endtime=10**50, species_selection=species_selection, track_propensities=self._IsTrackPropensities,rate_selection = rate_selection,last_timepoint=False,seed = self._IsSeed,quiet=quiet)      
        self.SSA.Execute(self.settings)     
        
        if self.settings.species_selection:   
            self.sim_species_tracked = [s_id for s_id in self.settings.species_selection]
        else:
            self.sim_species_tracked = copy.copy(self.SSA.species_names)           
                        
        (L_probability_mass, D_means, D_stds,D_moments) = Analysis.GetSpeciesDistributions(self.SSA.sim_output,self.sim_species_tracked)      
        m1 = [np.array(list(D_moments[s_id].values())) for s_id in self.sim_species_tracked]      
        IsContinue = True
        if not quiet:
            print('Info: {0:d} time steps simulated'.format(size))
        n=1      
        while IsContinue:          
            self.settings = SSASettings(x_matrix=self.SSA.X_matrixinit,timesteps=size*(n+1),starttime=self.SSA.sim_t,endtime=10**50, species_selection=species_selection, track_propensities=self._IsTrackPropensities,rate_selection = rate_selection,last_timepoint=False,seed = self._IsSeed,quiet=quiet)      
            self.SSA.Execute(self.settings)
            (L_probability_mass,D_means,D_stds,D_moments) = Analysis.GetDataDistributions(self.SSA.sim_output,self.sim_species_tracked)          
            m2 = [np.array(list(D_moments[s_id].values())) for s_id in self.sim_species_tracked] 
            max_total = 0
            for i in range(self.SSA.n_species): 
                max_s = abs(1-(m2[i]/m1[i])).max()
                if max_s > max_total:
                     max_total = max_s          
            m1 = copy.deepcopy(m2)  
            n+=1
            if not quiet:
                print('Info: {0:d} time steps simulated'.format(n*size))
            if max_total < error:
                IsContinue = False                  
        t2 = time.time()
        self.simulation_time = t2-t1
        if not quiet: 
            print("Info: Simulation time {0:1.5f}".format(self.simulation_time))
        self.FillDataStochsim()
        self._IsSimulationDone = True
        self.sim_trajectories_done = copy.copy(self.sim_trajectories)
        try:
            self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.sim_rates_tracked,self.plot.plotnum,quiet)
        except:
            self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.sim_rates_tracked,quiet=quiet)
            

    def _getSpecies2Plot(self,species2plot):
        """ *** For internal use only ***: this function determines the species for which we will plot something """
        if species2plot == True:
            species2plot = self.sim_species_tracked
        if isinstance(species2plot,str):
            species2plot = [species2plot]
        for s_id in species2plot:
            assert s_id in self.sim_species_tracked, "Species {0} is not in the model or species selection".format(s_id)
        return species2plot
        

    def _getRates2Plot(self,rates2plot):
        """ *** For internal use only ***: this function determines the reactions for which we will plot something """
        if rates2plot == True:
            rates2plot = self.sim_rates_tracked
        if isinstance(rates2plot,str):
            rates2plot = [rates2plot]
        for r_id in rates2plot:
            assert r_id in self.sim_rates_tracked, "Reaction {0} is not in the model or reaction selection".format(r_id)
        return rates2plot
        

    def GetWaitingtimes(self):
        """ Get for each reaction the waiting times """      
        assert self._IsSimulationDone, "First do a stochastic simulation (and do not use the Tau-Leaping method)"         
        assert not self._IsTauleaping, "Tau-Leaping method does not allow for calculation of waiting times"          
        assert not self._IsFixedIntervalMethod, "Fixed-interval output solvers do not allow for calculation of waiting times"
        assert not self._IsOnlyLastTimepoint, "Calculating waiting times is disabled when saving only the last time point"        
        
        for n in range(1,self.sim_trajectories_done+1): 
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            D_waitingtimes = Analysis.ObtainWaitingtimes(self.data_stochsim,self.SSA.rate_names) # hard coded for all reactions
            self.data_stochsim.setWaitingtimes(D_waitingtimes,self.SSA.rate_names)
            self.data_stochsim.setWaitingtimesMeans(self.data_stochsim.waiting_times,self.SSA.rate_names)        
            self.data_stochsim.setWaitingtimesStandardDeviations(self.data_stochsim.waiting_times,self.SSA.rate_names)
            if self.sim_trajectories_done > 1: # "store" the data, otherwise the added waiting times get lost again by opening via GetTrajectoryData
                self.DumpTrajectoryData(n)            
        

    def GetRegularGrid(self,n_samples=51):
        """
        The Gillespie method generates data at irregular time points. This function puts the data on a fixed regular time grid of which the user can specify the resolution (n_samples). 
        
        For each trajectory, we use the same grid. This has the following consequences for the two type of simulation modes:
        - time: the end time of each trajectory is identical
        - time steps: the end time of each trajectory is different. We select the minimal end time of each of these simulations and ignore the period afterwards.
        
        Input:
         - *n_samples* [default = 51] (integer)
        """
        assert self._IsSimulationDone, "First do a stochastic simulation"        
        assert not self._IsOnlyLastTimepoint, "Generating a regular grid is disabled when saving only the last time point"
        
        ### Determine the number of samples ###
        if isinstance(n_samples,int):
            pass
        elif type(n_samples) in [float,np.float64,np.float32]:
            print("*** WARNING ***: 'n_samples' must be an integer rather than a float; float {0} is rounded to {1:d}".format(n_samples,int(n_samples)))
            n_samples = int(n_samples)
        elif n_samples == True:
            n_samples = int(self.data_stochsim.simulation_endtime)    
        else:
            raise TypeError("Argument of GetRegularGrid() must be an integer")               
        self._n_samples = n_samples
        
        n_species = len(self.data_stochsim.species_labels)
        L_species = [[] for i in range(n_species)]
        if self._IsCellDivision:
           L_volumes = []  # hard coded one sort of volume 
        if self._IsTrackPropensities:
            n_rates = len(self.sim_rates_tracked)
            L_propensities = [[] for j in range(n_rates)]
            self.data_stochsim_grid.propensities_autocorrelations = [[] for j in range(n_rates)]
        
        if self.sim_mode == 'time':
            sample_timepoints = np.linspace(0,self.sim_end,n_samples)            
        else: 
            simulation_endtimes = []
            for n in range(1,self.sim_trajectories_done+1): 
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                simulation_endtimes.append(self.data_stochsim.simulation_endtime)                             
            sample_timepoints = np.linspace(0,min(simulation_endtimes),n_samples)
        
        for n in range(1,self.sim_trajectories_done+1):               
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)           
            
            sample_indices = self.data_stochsim.time.searchsorted(sample_timepoints, side = 'right') - 1 # Get point before
            sampled_species = self.data_stochsim.species[sample_indices]        
            if self._IsCellDivision:                        
                sampled_volume = self.data_stochsim.volume[sample_indices]             
                
            ### Put data in grid files ### 
            for i in range(n_species):
                L_species[i].append(sampled_species[:,i])             
             
            if self._IsTrackPropensities:
                sample_propensities = self.data_stochsim.propensities[sample_indices]   
                for j in range(n_rates):
                    L_propensities[j].append(sample_propensities[:,j])     
            if self._IsCellDivision:
                L_volumes.append(sampled_volume)
                    
        self.data_stochsim_grid.setTime(sample_timepoints)                        
        self.data_stochsim_grid.setSpecies(L_species,self.sim_species_tracked)      
        (self.data_stochsim_grid.species_means,self.data_stochsim_grid.species_standard_deviations) = Analysis.GetAverageResults(self.data_stochsim_grid.species)
        if self._IsTrackPropensities:
            self.data_stochsim_grid.setPropensities(L_propensities,self.sim_rates_tracked)              
            (self.data_stochsim_grid.propensities_means,self.data_stochsim_grid.propensities_standard_deviations) = Analysis.GetAverageResults(self.data_stochsim_grid.propensities)
        if self._IsCellDivision:
            self.data_stochsim_grid.setVolume(L_volumes)
            (self.data_stochsim_grid.volume_means,self.data_stochsim_grid.volume_standard_deviations) = Analysis.GetAverageResults(self.data_stochsim_grid.volume)
        self.HAS_AVERAGE = True        
            

    def GetAverageSpeciesDistributions(self):
        """ Get average species distributions """      
        assert self._IsSimulationDone, "First do a stochastic simulation."         
        assert not self._IsOnlyLastTimepoint, "Calculating average species distributions is disabled when saving only the last time point"
        
        D_distributions = {s_id:{}  for s_id in self.sim_species_tracked}
        L_distributions_means = []
        L_distributions_standard_deviations = []
        for n in range(1,self.sim_trajectories_done+1): 
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            for i in range(len(self.sim_species_tracked)):
                s_id = self.sim_species_tracked[i]
                for m,s_amount in enumerate(self.data_stochsim.species_distributions[i][0]):
                    if not s_amount in list(D_distributions[s_id]):
                        D_distributions[s_id][s_amount] = []
                    D_distributions[s_id][s_amount].append(self.data_stochsim.species_distributions[i][1][m])
            
        for s_id in self.sim_species_tracked:
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
        self.data_stochsim_grid.setSpeciesDistributionAverage(L_distributions_means,L_distributions_standard_deviations)          
        

    def GetAveragePropensitiesDistributions(self):
        """ Get average propensities distributions """
        assert (self._IsSimulationDone and self._IsTrackPropensities), "First do a stochastic simulation (use the IsTrackPropensities flag in DoStochSim)."
        assert not self._IsOnlyLastTimepoint, "Calculating average propensities distributions is disabled when saving only the last time point"
        
        D_distributions = {r_id:{}  for r_id in self.sim_rates_tracked}
        L_distributions_means = []
        L_distributions_standard_deviations = []
        for n in range(1,self.sim_trajectories_done+1): 
            if self.sim_trajectories_done > 1:
                self.GetTrajectoryData(n)
            for j in range(len(self.sim_rates_tracked)):
                r_id = self.sim_rates_tracked[j]                
                for m,r_prop in enumerate(self.data_stochsim.propensities_distributions[j][0]):
                    if not r_prop in list(D_distributions[r_id]): 
                        D_distributions[r_id][r_prop] = []
                    D_distributions[r_id][r_prop].append(self.data_stochsim.propensities_distributions[j][1][m])
            
        for r_id in self.sim_rates_tracked:
            L_propensities = list(D_distributions[r_id])  # for a given species 
            L_means = []
            L_stds = []   
            for r_prop in L_propensities:
                while len(D_distributions[r_id][r_prop]) < (n-1):
                    D_distributions[r_id][r_prop].append(0)
                L_means.append(np.mean(D_distributions[r_id][r_prop]))
                L_stds.append(np.std(D_distributions[r_id][r_prop]))    
            L_distributions_means.append([L_propensities,L_means])
            L_distributions_standard_deviations.append([L_propensities,L_stds])
        self.data_stochsim_grid.setPropensitiesDistributionAverage(L_distributions_means,L_distributions_standard_deviations)            


    def GetSpeciesAutocorrelations(self,species2calc=True,n_samples=51):
        """
        Input:
         - *species2calc* [default = True] as a list ['S1','S2']
         - *n_samples* (integer)
        """
        if (not self.HAS_AVERAGE) or n_samples != 51:
            self.GetRegularGrid(n_samples)        
        
        species2calc = self._getSpecies2Plot(species2calc)                             
        L_species_autocorrelations = [[] for i in range(len(self.sim_species_tracked))]        
        for n in range(self.sim_trajectories_done):               
            for s_id in species2calc:                  
                s_index = self.sim_species_tracked.index(s_id)
                L_species_autocorrelations[s_index].append(Analysis.Autocorrelation(self.data_stochsim_grid.species[s_index][n]))  
            
        self.data_stochsim_grid.setSpeciesAutocorrelations(L_species_autocorrelations)
        
        (self.data_stochsim_grid.species_autocorrelations_means,
        self.data_stochsim_grid.species_autocorrelations_standard_deviations) = Analysis.GetAverageResults(L_species_autocorrelations)


    def GetSpeciesAutocovariances(self,species2calc=True,n_samples=51):
        """
        Input:
         - *species2calc* [default = True] as a list ['S1','S2']
         - *n_samples* (integer)
        """
        if (not self.HAS_AVERAGE) or n_samples != 51:
            self.GetRegularGrid(n_samples)
        
        species2calc = self._getSpecies2Plot(species2calc)                             
        L_species_autocovariances = [[] for i in range(len(self.sim_species_tracked))]       
        for n in range(self.sim_trajectories_done):                  
            for s_id in species2calc:
                s_index = self.sim_species_tracked.index(s_id)  
                L_species_autocovariances[s_index].append(Analysis.AutoCov(self.data_stochsim_grid.species[s_index][n]))  
                  
        self.data_stochsim_grid.setSpeciesAutocovariances(L_species_autocovariances)
        (self.data_stochsim_grid.species_autocovariances_means,
        self.data_stochsim_grid.species_autocovariances_standard_deviations) = Analysis.GetAverageResults(L_species_autocovariances)
            

    def GetPropensitiesAutocorrelations(self,rates=True,n_samples=51):
        """
        Input:
         - *rates* [default = True] as a list ['R1','R2']     
         - *n_samples* (integer)
        """
        assert self._IsTrackPropensities, "First do a stochastic simulation with tracking propensities (IsTrackPropensities=True in DoStochSim)"       
        
        if (not self.HAS_AVERAGE) and n_samples != 51:
            self.GetRegularGrid(n_samples)          
        
        rates = self._getRates2Plot(rates)
        L_Propensities_autocorrelations = [[] for j in range(len(self.sim_rates_tracked))]        
        for n in range(self.sim_trajectories_done):               
            for r_id in rates:
                r_index = self.sim_rates_tracked.index(r_id)   
                L_Propensities_autocorrelations[r_index].append(Analysis.Autocorrelation(self.data_stochsim_grid.propensities[r_index][n]))
                    
        self.data_stochsim_grid.setPropensitiesAutocorrelations(L_Propensities_autocorrelations)
        self.data_stochsim_grid.propensities_autocorrelations_means, self.data_stochsim_grid.propensities_autocorrelations_standard_deviations = Analysis.GetAverageResults(self.data_stochsim_grid.propensities_autocorrelations)


    def GetPropensitiesAutocovariances(self,rates=True,n_samples=51): 
        """
        Input:
         - *rates* [default = True] as a list ['R1','R2']     
         - *n_samples* (integer)
        """
        assert self._IsTrackPropensities, "First do a stochastic simulation with tracking propensities (IsTrackPropensities=True in DoStochSim)"
        
        if (not self.HAS_AVERAGE) and n_samples != 51:
            self.GetRegularGrid(n_samples)
            
        rates = self._getRates2Plot(rates)
        L_Propensities_autocovariances = [[] for j in range(len(self.sim_rates_tracked))]              
        for n in range(self.sim_trajectories_done):                
            for r_id in rates:
                r_index = self.sim_rates_tracked.index(r_id)   
                L_Propensities_autocovariances[r_index].append(Analysis.AutoCov(self.data_stochsim_grid.propensities[r_index][n]))
              
        self.data_stochsim_grid.setPropensitiesAutocovariances(L_Propensities_autocovariances)
        self.data_stochsim_grid.propensities_autocovariances_means, self.data_stochsim_grid.propensities_autocovariances_standard_deviations = Analysis.GetAverageResults(self.data_stochsim_grid.propensities_autocovariances)                      
        

    def Export2File(self,analysis='timeseries',datatype='species', IsAverage = False, directory=None,quiet=False):
        """
        Write data to a text document     
    
        Input:
         - *analysis* [default = 'timeseries'] (string) options: timeseries, distribution, mean, std, autocorrelation, autocovariance
         - *datatype*  [default = 'species'] (string) options: species, propensities, waitingtimes
         - *IsAverage* [default = False] (boolean)   
         - *directory* [default = None] (string)
         - *quiet* [default = False] (boolean)
        """
        if self._IsQuiet:
            quiet = True        
        if directory == None:
            if not IsAverage:
                directory = os.path.join(self.output_dir,"{0:s}_{1:s}_{2:s}".format(self.model_file,datatype,analysis))
            else:
                directory = os.path.join(self.output_dir,"{0:s}_{1:s}_{2:s}_{3:s}".format(self.model_file,"average",datatype,analysis))
        else:
            if not os.path.exists(directory):
                os.makedirs(directory)
            if not IsAverage:
                directory = os.path.join(directory,"{0:s}_{1:s}_{2:s}".format(self.model_file,datatype,analysis))  
            else:
                directory = os.path.join(directory,"{0:s}_{1:s}_{2:s}_{3:s}".format(self.model_file,"average",datatype,analysis))
          
        if (datatype.lower() == 'species') and (analysis.lower() == 'timeseries') and (not IsAverage):
            assert self._IsSimulationDone, "First do a stochastic simulation"          
            for n in range(1,self.sim_trajectories_done+1): 
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)	# Dir/Filename
                file_out = open(file_path,'w')
                file_out.write('# reactions: {0:d}\n'.format(len(self.SSA.rate_names)))
                file_out.write('Time')
                for s_id in self.data_stochsim.species_labels:
                    file_out.write('\t{0:s}'.format(s_id))
                
                if not self._IsTauleaping:
                    file_out.write('\tFired Reaction')              
                file_out.write('\n')                                                       
                for m,timepoint in enumerate(self.data_stochsim.getSpecies()):
                    slist = [str(value) for value in timepoint]
                    line = "\t".join(slist) 
                    if not self._IsTauleaping:
                        line += "\t {0}".format(self.data_stochsim.fired_reactions[m])
                    line += '\n'
                    file_out.write(line)                
                file_out.close()  
                if not quiet:              
                    print("Info: Species time series output is successfully saved at: {0:s}".format(file_path) )
                  
        elif (datatype.lower() == 'species') and (analysis.lower() in ['distributions','distribution'] ) and (not IsAverage):
            assert self._IsSimulationDone, "First do a stochastic simulation"          
            for n in range(1,self.sim_trajectories_done+1):  
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')  
                for L_species_dist in self.data_stochsim.species_distributions:
                    file_out.write("Copy number\tPMF\n")
                    for m in range(len(L_species_dist[0])):
                        file_out.write("{0}\t{1}\n".format(L_species_dist[0][m],L_species_dist[1][m]) ) 
                file_out.close()
                if not quiet:
                    print("Info: Species distributions output is successfully saved at: {0:s}".format(file_path) )
                  
        elif  (datatype.lower() == 'species') and (analysis.lower() in ['autocorrelation','autocorrelations']) and (not IsAverage):
            if not self.data_stochsim_grid.HAS_SPECIES_AUTOCORRELATIONS:
                print("*** WARNING ***: Autocorrelations are not yet calculated. StochPy automatically calculates autocorrelations with pre-defined settings. You can use GetSpeciesAutocorrelations(species2calc=True,n_samples=51)")
                self.GetSpeciesAutocorrelations()          
             
            for n in range(1,self.sim_trajectories_done+1):
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                  
                Arr_autocorrelations = self.data_stochsim_grid.species_autocorrelations[:,n-1]
                file_out.write('Lag (tau)')
                for s_id in self.sim_species_tracked:
                    file_out.write('\tAutocorrelation ({0:s})'.format(s_id))
                file_out.write('\n')                 
                for x,t in enumerate(self.data_stochsim_grid.time):               
                    file_out.write(str(t)) # lag                      
                    for acor in Arr_autocorrelations[:,x]:
                        file_out.write('\t{0:f}'.format(acor))
                    file_out.write('\n')                                
                               
                file_out.close()
                if not quiet:
                    print("Info: Species autocorrelation output is successfully saved at: {0:s}".format(file_path) )

        elif  (datatype.lower() == 'species') and (analysis.lower() in ['autocovariance','autocovariances']) and (not IsAverage):
            if not self.data_stochsim_grid.HAS_SPECIES_AUTOCOVARIANCES:
                print("*** WARNING ***: Autocovariances are not yet calculated. StochPy automatically calculates autocovariances with pre-defined settings. You can use GetSpeciesAutocovariances(species2calc=True,n_samples=51)")
                self.GetSpeciesAutocovariances()          
            for n in range(1,self.sim_trajectories_done+1):
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                  
                Arr_autocovariances = self.data_stochsim_grid.species_autocovariances[:,n-1]
                file_out.write('Lag (tau)')
                for s_id in self.sim_species_tracked:
                    file_out.write('\tAutocovariance ({0:s})'.format(s_id))
                file_out.write('\n')                 
                for x,t in enumerate(self.data_stochsim_grid.time):               
                    file_out.write(str(t)) # lag                      
                    for acov in Arr_autocovariances[:,x]:
                        file_out.write('\t{0:f}'.format(acov))
                    file_out.write('\n')                    
                file_out.close()
                if not quiet:
                    print("Info: Species autocovariance output is successfully saved at: {0:s}".format(file_path) )
                    
        elif (datatype.lower() == 'species') and (analysis.lower() == 'mean') and (not IsAverage):
            assert self._IsSimulationDone, "First do a stochastic simulation"          
            for n in range(1,self.sim_trajectories_done+1):
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                      
                file_out.write("Species\tMean\n")                
                for s_id in self.sim_species_tracked: 
                    file_out.write("{0:s}\t{1:f}\n".format(s_id,self.data_stochsim.species_means[s_id])) 
                file_out.close()
                if not quiet:
                    print("Info: Species means output is successfully saved at: {0:s}".format(file_path) )
                    
        elif (datatype.lower() == 'species') and (analysis.lower() == 'std') and (not IsAverage):
            assert self._IsSimulationDone, "First do a stochastic simulation"          
            for n in range(1,self.sim_trajectories_done+1):
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                      
                file_out.write("Species\tStandard deviation\n")                
                for s_id in self.sim_species_tracked: 
                    file_out.write("{0:s}\t{1:f}\n".format(s_id,self.data_stochsim.species_standard_deviations[s_id]))               
                file_out.close()
                if not quiet:
                    print("Info: Species standard deviations output is successfully saved at: {0:s}".format(file_path) )
                    
        elif (datatype.lower() == 'propensities') and  (analysis.lower() == 'timeseries') and (not IsAverage):       
            assert (self._IsTrackPropensities and self._IsSimulationDone), "First do a stochastic simulation with tracking propensities (IsTrackPropensities=True in DoStochSim)"
            for n in range(1,self.sim_trajectories_done+1):
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')  
                file_out.write('Time')
                for r_id in self.sim_rates_tracked:
                    file_out.write('\t{0:s}'.format(r_id) )
                file_out.write('\n')      
                for timepoint in self.data_stochsim.getPropensities(): 
                    slist = [str(value) for value in timepoint]
                    line = "\t".join(slist) 
                    line += '\n'
                    file_out.write(line)                
                file_out.close()
                if not quiet:
                    print("Info: Propensities time series output is successfully saved at: {0:s}".format(file_path) )
                    
        elif (datatype.lower() == 'propensities') and (analysis.lower() in ['distributions','distribution'] ) and (not IsAverage):
            assert (self._IsTrackPropensities and self._IsSimulationDone), "First do a stochastic simulation with tracking propensities (IsTrackPropensities=True in DoStochSim)"
            for n in range(1,self.sim_trajectories_done+1):
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                
                for j,L_prop_dist in enumerate(self.data_stochsim.propensities_distributions): # TODO: stuff
                    file_out.write("Propensity ({0:s})\tPMF\n".format(self.sim_rates_tracked[j]) )
                    for k in range(len(L_prop_dist[0])):
                        file_out.write("{0}\t{1}\n".format(L_prop_dist[0][k],L_prop_dist[1][k]))                
                file_out.close()
                if not quiet:
                    print("Info: Species distributions output is successfully saved at: {0:s}".format(file_path) )
                
        elif (datatype.lower() == 'propensities') and (analysis.lower() in ['autocorrelation','autocorrelations'] ) and (not IsAverage): # TODO: test
            if not self.data_stochsim_grid.HAS_PROPENSITIES_AUTOCORRELATIONS:
                print("*** WARNING ***: Autocorrelations are not yet calculated. StochPy automatically calculates autocorrelations with pre-defined settings. You can use GetPropensitiesAutocorrelations(rates=True,n_samples=51)")
                self.GetPropensitiesAutocorrelations(n_samples=self.data_stochsim.simulation_endtime)          
            for n in range(1,self.sim_trajectories_done+1):
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')  
                Arr_autocorrelations = self.data_stochsim_grid.propensities_autocorrelations[:,n-1]
                file_out.write('Lag (tau)')
                for r_id in self.sim_rates_tracked:
                    file_out.write('\tAutocorrelation ({0:s})'.format(r_id))
                file_out.write('\n')
                for x,t in enumerate(self.data_stochsim_grid.time):               
                    file_out.write(str(t)) # lag                      
                    for acor in Arr_autocorrelations[:,x]:
                        file_out.write('\t{0:f}'.format(acor))
                    file_out.write('\n')              
            
                file_out.close()
                if not quiet:
                    print("Info: Propensities autocorrelation output is successfully saved at: {0:s}".format(file_path) ) 

        elif (datatype.lower() == 'propensities') and (analysis.lower() in ['autocovariance','autocovariances']) and (not IsAverage):
            if not self.data_stochsim_grid.HAS_PROPENSITIES_AUTOCOVARIANCES:
                print("*** WARNING ***: Autocovariances are not yet calculated. StochPy automatically calculates autocovariances with pre-defined settings. You can use GetPropensitiesAutocovariances(rates=True,n_samples=51)")
                self.GetPropensitiesAutocovariances(rates = rates2plot,n_samples=self.data_stochsim.simulation_endtime)
            
            for n in range(1,self.sim_trajectories_done+1): 
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')  
                Arr_autocovariances = self.data_stochsim_grid.propensities_autocovariances[:,n-1]
                file_out.write('Lag (tau)')
                for r_id in self.sim_rates_tracked:
                    file_out.write('\tAutocovariance ({0:s})'.format(r_id) )          
                file_out.write('\n')
                for x,t in enumerate(self.data_stochsim_grid.time):               
                    file_out.write(str(t)) # lag                      
                    for acov in Arr_autocovariances[:,x]:
                        file_out.write('\t{0:f}'.format(acov))
                    file_out.write('\n')   
                                       
                file_out.close()
                if not quiet:
                    print("Info: Propensities autocovariance output is successfully saved at: {0:s}".format(file_path) )
                  
        elif (datatype.lower() == 'propensities') and  (analysis.lower() == 'mean') and (not IsAverage):        
            assert (self._IsTrackPropensities and self._IsSimulationDone), "First do a stochastic simulation with tracking propensities (IsTrackPropensities=True in DoStochSim)"
            for n in range(1,self.sim_trajectories_done+1):                
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                      
                file_out.write("Reaction\tMean\n")                
                for r_id in self.sim_rates_tracked: 
                    file_out.write("{0}\t{1}\n".format(r_id,self.data_stochsim.propensities_means[r_id]))                                    
                file_out.close()
                if not quiet:
                    print("Info: Propensities means output is successfully saved at: {0:s}".format(file_path) )
                                    
        elif (datatype.lower() == 'propensities') and (analysis.lower() == 'std') and (not IsAverage):
            assert (self._IsTrackPropensities and self._IsSimulationDone), "First do a stochastic simulation with tracking propensities (IsTrackPropensities=True in DoStochSim)"
            for n in range(1,self.sim_trajectories_done+1):            
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                      
                file_out.write("Reaction\tStandard deviation\n")                
                for r_id in self.sim_rates_tracked: 
                    file_out.write("{0}\t{1}\n".format(r_id,self.data_stochsim.propensities_standard_deviations[r_id]))                
                file_out.close()
                if not quiet:
                    print("Info: Propensities standard deviations output is successfully saved at: {0:s}".format(file_path) )
                  
        elif (datatype.lower() == 'waitingtimes') and (analysis.lower() in ['distributions','distribution'] ) and (not IsAverage):          
            if (not self.data_stochsim.HAS_WAITINGTIMES):
                self.GetWaitingtimes()
            for n in range(1,self.sim_trajectories_done+1):   
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')
                for r_id in sorted(self.data_stochsim.waiting_times):
                    file_out.write("Waitingtimes\t{0:s}\n".format(r_id) )
                    L_waiting_times_r = self.data_stochsim.waiting_times[r_id]
                    for time in L_waiting_times_r:
                        file_out.write("{0:f}\n".format(time) )                
                file_out.close()
                if not quiet:
                    print("Info: Waitingtimes distributions output is successfully saved at: {0:s}".format(file_path) )
                  
        elif (datatype.lower() == 'waitingtimes') and (analysis.lower() == 'mean') and (not IsAverage):
            if (not self.data_stochsim.HAS_WAITINGTIMES):
                self.GetWaitingtimes()          
            for n in range(1,self.sim_trajectories_done+1):              
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                      
                file_out.write("Reaction\tMean\n")                
                for j,r_id in enumerate(self.SSA.rate_names): 
                    file_out.write("{0}\t{1}\n".format(r_id,self.data_stochsim.waiting_times_means[j]))                 
                file_out.close()
                if not quiet:
                    print("Info: Waitingtimes means output is successfully saved at: {0:s}".format(file_path) )
                  
        elif (datatype.lower() == 'waitingtimes') and (analysis.lower() == 'std') and (not IsAverage):
            if (not self.data_stochsim.HAS_WAITINGTIMES):
                self.GetWaitingtimes()          
            for n in range(1,self.sim_trajectories_done+1):            
                if self.sim_trajectories_done > 1:
                    self.GetTrajectoryData(n)
                file_path = "{0:s}{1:d}.txt".format(directory,n)
                file_out = open(file_path,'w')                      
                file_out.write("Reaction\tStandard deviation\n")                
                for j,r_id in enumerate(self.SSA.rate_names): 
                    file_out.write("{0}\t{1}\n".format(r_id,self.data_stochsim.waiting_times_standard_deviations[j]))                                    
                file_out.close()
                if not quiet:
                    print("Info: Waitingtimes means output is successfully saved at: {0:s}".format(file_path) )  
                   
        elif (datatype.lower() == 'species') and (analysis.lower() == 'timeseries') and (IsAverage):
            if not self.HAS_AVERAGE:
                print("*** WARNING ***: No regular grid is created yet. Use GetRegularGrid(n_samples) if averaged results are unsatisfactory")
                self.GetRegularGrid()          
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')
            file_out.write("t")
            for s_id in self.sim_species_tracked:
                file_out.write("\t{0:s} (Mean)\t{0:s} (STD)".format(s_id) )
            file_out.write("\n")            
            for x,t in enumerate(self.data_stochsim_grid.time): 
                file_out.write( str(t) )  
                for i in range(len(self.sim_species_tracked)): 
                    file_out.write("\t{0:f}\t{1:f}".format(self.data_stochsim_grid.species_means[x,i],self.data_stochsim_grid.species_standard_deviations[x,i]) )
                file_out.write("\n")
            if not quiet:    
                print("Info: Averaged species time series output is successfully saved at: {0:s}".format(file_path) )
            
        elif (datatype.lower() == 'propensities') and (analysis.lower() == 'timeseries') and (IsAverage):
            assert self._IsTrackPropensities, "First do a stochastic simulation with tracking propensities (IsTrackPropensities=True in DoStochSim)"
            if not self.HAS_AVERAGE and self._IsTrackPropensities:
                print("*** WARNING ***: No regular grid is created yet. Use GetRegularGrid(n_samples) if averaged results are unsatisfactory")
                self.GetRegularGrid()          
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')
            file_out.write("t")
            for r_id in self.sim_rates_tracked:
                file_out.write("\t{0:s} (Mean)\t{0:s} (STD)".format(r_id) )
            file_out.write("\n")                     
            for x,t in enumerate(self.data_stochsim_grid.time):
                file_out.write( str(t) )
                for j in range(len(self.sim_rates_tracked)): 
                    file_out.write("\t{0:f}\t{1:f}".format(self.data_stochsim_grid.propensities_means[x,j],self.data_stochsim_grid.propensities_standard_deviations[x,j]) )
                file_out.write("\n") 
            if not quiet:                   
                print("Info: Averaged propensities time series output is successfully saved at: {0:s}".format(file_path) )
            
        elif (datatype.lower() == 'species') and (analysis.lower() in ['distributions','distribution'] ) and (IsAverage):
            if not self.data_stochsim_grid.HAS_AVERAGE_SPECIES_DISTRIBUTIONS:              
                self.GetAverageSpeciesDistributions()          
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')             
            for i,s_id in enumerate(self.sim_species_tracked):
                file_out.write("Copy number\t{0:s} (Mean)\t{0:s} (STD)\n".format(s_id) )   
                for m in range(len(self.data_stochsim_grid.species_distributions_means[i][0])):
                    s_amount = self.data_stochsim_grid.species_distributions_means[i][0][m] 
                    s_probability_mean = self.data_stochsim_grid.species_distributions_means[i][1][m]
                    s_probability_std = self.data_stochsim_grid.species_distributions_standard_deviations[i][1][m]
                    file_out.write("{0:0.0f}\t{1:f}\t{2:f}\n".format(s_amount,s_probability_mean,s_probability_std) )
                file_out.write("\n")   
            if not quiet:                 
                print("Info: Averaged species distributions output is successfully saved at: {0:s}".format(file_path) )
            
        elif (datatype.lower() == 'propensities') and (analysis.lower() in ['distributions','distribution'] ) and (IsAverage):
            if not self.data_stochsim_grid.HAS_AVERAGE_PROPENSITIES_DISTRIBUTIONS:
                self.GetAveragePropensitiesDistributions()
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')            
            for j,r_id in enumerate(self.sim_rates_tracked):
                file_out.write("Propensity\t{0:s} (Mean)\t{0:s} (STD)\n".format(r_id) )   
                for m in range(len(self.data_stochsim_grid.propensities_distributions_means[j][0])):
                    r_prop = self.data_stochsim_grid.propensities_distributions_means[j][0][m]
                    r_probability_mean = self.data_stochsim_grid.propensities_distributions_means[j][1][m]
                    r_probability_std = self.data_stochsim_grid.propensities_distributions_standard_deviations[j][1][m] 
                    file_out.write("{0}\t{1:f}\t{2:f}\n".format(r_prop,r_probability_mean,r_probability_std))
                file_out.write("\n")
            if not quiet:               
                print("Info: Averaged propensities distributions output is successfully saved at: {0:s}".format(file_path) )
            
        elif (datatype.lower() == 'species') and (analysis.lower() in ['autocorrelation','autocorrelations'] ) and (IsAverage):         
            if not self.data_stochsim_grid.HAS_SPECIES_AUTOCORRELATIONS:
                print("*** WARNING ***: Autocorrelations are not yet calculated. StochPy automatically calculates autocorrelations with pre-defined settings. You can use GetSpeciesAutocorrelations(species2calc=True,n_samples=51)")
                self.GetSpeciesAutocorrelations()
            
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')                        
            file_out.write("t")
            for s_id in self.sim_species_tracked:
                file_out.write("\t{0:s} (Mean)\t{0:s} (STD)".format(s_id) )
            file_out.write("\n")                  

            for x,t in enumerate(self.data_stochsim_grid.time):                
                file_out.write(str(t))
                for i in range(len(self.sim_species_tracked)): 
                    acor_mean = self.data_stochsim_grid.species_autocorrelations_means[x,i]
                    acor_std = self.data_stochsim_grid.species_autocorrelations_standard_deviations[x,i]  
                    file_out.write("\t{0:f}\t{1:f}".format(acor_mean,acor_std) )
                file_out.write("\n")                               
                
            if not quiet:                  
                print("Info: Averaged species autocorrelations output is successfully saved at: {0:s}".format(file_path) )
            
        elif (datatype.lower() == 'species') and (analysis.lower() in ['autocovariance','autocovariances']) and (IsAverage):          
            if not self.data_stochsim_grid.HAS_SPECIES_AUTOCOVARIANCES:
                print("*** WARNING ***: Autocovariances are not yet calculated. StochPy automatically calculates autocovariances with pre-defined settings. You can use GetSpeciesAutocovariances(species2calc=True,n_samples=51)")
                self.GetSpeciesAutocovariances()
            
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')                        
            file_out.write("t")
            for s_id in self.sim_species_tracked:
                file_out.write("\t{0:s} (Mean)\t{0:s} (STD)".format(s_id) )
            file_out.write("\n")                  

            for x,t in enumerate(self.data_stochsim_grid.time):                
                file_out.write(str(t))
                for i in range(len(self.sim_species_tracked)): 
                    acor_mean = self.data_stochsim_grid.species_autocovariances_means[x,i]
                    acor_std = self.data_stochsim_grid.species_autocovariances_standard_deviations[x,i]  
                    file_out.write("\t{0:f}\t{1:f}".format(acor_mean,acor_std) )
                file_out.write("\n")          
            if not quiet:              
                print("Info: Averaged species autocovariances output is successfully saved at: {0:s}".format(file_path) )
            
        elif (datatype.lower() == 'propensities') and (analysis.lower() in ['autocorrelation','autocorrelations']) and (IsAverage):    # TODO: test       
            if not self.data_stochsim_grid.HAS_PROPENSITIES_AUTOCORRELATIONS:
                print("*** WARNING ***: Autocorrelations are not yet calculated. StochPy automatically calculates autocorrelations with pre-defined settings. You can use GetPropensitiesAutocorrelations(rates=True,n_samples=51)")
                self.GetPropensitiesAutocorrelations()          
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')  
            
            file_out.write("t")
            for r_id in self.sim_rates_tracked:
                file_out.write("\t{0:s} (Mean)\t{0:s} (STD)".format(r_id) )
            file_out.write("\n")

            for x,t in enumerate(self.data_stochsim_grid.time):                
                file_out.write(str(t))
                for j in range(len(self.sim_rates_tracked)): 
                    acor_mean = self.data_stochsim_grid.propensities_autocorrelations_means[x,j]
                    acor_std = self.data_stochsim_grid.propensities_autocorrelations_standard_deviations[x,j]  
                    file_out.write("\t{0:f}\t{1:f}".format(acor_mean,acor_std) )
                file_out.write("\n")              
            if not quiet:               
                print("Info: Averaged propensities autocorrelations output is successfully saved at: {0:s}".format(file_path) )

        elif (datatype.lower() == 'propensities') and (analysis.lower() in ['autocovariance','autocovariances']) and (IsAverage):       # TODO: test   
            if not self.data_stochsim_grid.HAS_PROPENSITIES_AUTOCOVARIANCES:
                print("*** WARNING ***: Autocovariances are not yet calculated. StochPy automatically calculates autocovariances with pre-defined settings. You can use GetPropensitiesAutocovariances(rates=True,n_samples=51)")
                self.GetPropensitiesAutocovariances()
            
            file_path = '{0:s}.txt'.format(directory)
            file_out = open(file_path,'w')                        
            file_out.write("t")
            for r_id in self.sim_rates_tracked:
                file_out.write("\t{0:s} (Mean)\t{0:s} (STD)".format(r_id) )
            file_out.write("\n")                  

            for x,t in enumerate(self.data_stochsim_grid.time):                
                file_out.write(str(t))
                for i in range(len(self.sim_species_tracked)): 
                    acor_mean = self.data_stochsim_grid.propensities_autocovariances_means[x,i]
                    acor_std = self.data_stochsim_grid.propensities_autocovariances_standard_deviations[x,i]  
                    file_out.write("\t{0:f}\t{1:f}".format(acor_mean,acor_std) )
                file_out.write("\n")   
            if not quiet:              
                print("Info: Averaged propensities autocovariances output is successfully saved at: {0:s}".format(file_path) )

        else:
            raise UserWarning("No valid option specified. Nothing is exported. See help function (help(Export2File))")
            

    def Import2StochPy(self,filename,filedir,delimiter='\t'):
        """
        Can import time series data with the following format:
        
        Time S1 S2 S3 Fired Reactions
        0    1  0  1  nan
        1.5  0  0  2  1
        etc.
        
        or 
        
        Time S1 S2 S3
        0    1  0  1
        1.5  0  0  2
        etc.      
        
        In the future, this function will probably support more data types. We currently accept the default output of the Gillespie algorithm from which other data types can be derived.
        
        Input: 
         - *filename* (string)
         - *filedir* (string)
         - *delimiter* (string)
        """
        print("*** WARNING ***: In construction - This function currently accepts species time series data only!")        
        try:
            filepath = os.path.join(filedir,filename)
            file_in = open(filepath,'r') 
        except IOError:
            print("File path {0:s} does not exist".format(filepath) )
            sys.exit()
        
        L_sim_output = []
        L_data = file_in.readlines()
          
        nreactions = int(L_data[0].strip().split(':')[1])
        self.SSA.rate_names = []
        for r in range(1,nreactions+1):
            self.SSA.rate_names.append('R{0:d}'.format(r) )          

        print("Info: Reaction identifiers are unknown, StochPy created the following identifiers automatically:\t{0}".format(self.SSA.rate_names) )
        print("Info: Update 'smod.SSA.rate_names' if other identifiers are desired")                
        
        L_header = L_data[1].strip().split(delimiter)
        L_species_names = L_header[1:]      
        self._IsTauleaping = False
        if L_header[-1] != 'Fired Reaction': # no fired reactions stored                    
            self._IsTauleaping = True
        else:
            L_species_names.pop(-1)         
             
        for dat in L_data[2:]:
            fdat = [float(x) for x in dat.strip().split(delimiter)] # parse data and convert to float (is for some reason faster than integers ...)
            if not self._IsTauleaping and not math.isnan(fdat[-1]):
               fdat[-1] = int(fdat[-1])                  
            L_sim_output.append(fdat)      
        self.SSA.sim_output = L_sim_output
        self.SSA.species_names = L_species_names
        self.SSA.timestep = len(L_sim_output)
        self.SSA.sim_t = L_sim_output[-1][0]

        self._current_trajectory = 1 # HARD CODED
        self.sim_trajectories_done = 1      
        self._IsTrackPropensities =  False
        self.data_stochsim = IntegrationStochasticDataObj()           
        self.FillDataStochsim(IsImport=True)     
        try:
            self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.SSA.rate_names,self.plot.plotnum,quiet)
        except:
            self.plot = Analysis.DoPlotting(self.data_stochsim.species_labels,self.SSA.rate_names,quiet=quiet)
        self._IsSimulationDone = True        
            

    def ShowSpecies(self):
        """ Print the species of the model """
        print(self.SSA.species_names)


    def ShowOverview(self):
        """ Print an overview of the current settings """
        print("Current Model:\t{0:s}".format(self.model_file) )
        if self.sim_mode == "steps": 
            print("Number of time steps:\t{0:d}".format(self.sim_end) )
        elif self.sim_mode == "time":
            print("Simulation end time:\t{0:1.2f}".format(self.sim_end) )
        print("Current Algorithm:\t{0:s}".format(self.sim_method) )
        print("Number of trajectories:\t{0:d}".format(self.sim_trajectories) )
        if self._IsTrackPropensities:
             print("Propensities are tracked")
        else:
             print("Propensities are not tracked")


    def DeleteTempfiles(self):
        """ Deletes all .dat files """
        for name in os.listdir(self.temp_dir):
            try:
                os.remove(os.path.join(self.temp_dir,name))
            except:
                dir2delete = os.path.join(self.temp_dir,name)
                shutil.rmtree(dir2delete, ignore_errors=True)


    def DoTestsuite(self,epsilon_ = 0.01,sim_trajectories=1000):
        """
        DoTestsuite(epsilon_ = 0.01,sim_trajectories=1000)
        
        Do "sim_trajectories" simulations until t=50 and print the interpolated result for t = 0,1,2,...,50
        
        Input:
         - *epsilon_* [default = 0.01]: useful for tau-Leaping simulations (float)
         - *sim_trajectories* [default = 1000]
        """      
        self.DoStochSim(mode='time',end=50,epsilon = epsilon_,trajectories = sim_trajectories)
        self.GetRegularGrid(n_samples = 51)
        self.PrintAverageSpeciesTimeSeries()        


    def FillDataStochsim(self,IsImport=False):
        """
        Put all simulation data in the data object data_stochsim
        
        Input:
         - *IsImport* [default = False] (boolean)
        """ 
        if not IsImport and self.settings.species_selection:   
            self.sim_species_tracked = [s_id for s_id in self.settings.species_selection]            
        else:
            self.sim_species_tracked = copy.copy(self.SSA.species_names)

        if not IsImport and self.settings.rate_selection:
            self.sim_rates_tracked = [r_id for r_id in self.settings.rate_selection]
        else:
            self.sim_rates_tracked = copy.copy(self.SSA.rate_names)             
        
        (L_probability_mass,D_means,D_stds,D_moments) = Analysis.GetSpeciesDistributions(self.SSA.sim_output,self.sim_species_tracked)    
        sim_dat = np.array(self.SSA.sim_output)   
        self.data_stochsim.setTime(sim_dat[:,0])
        self.data_stochsim.setSpeciesDistributions(L_probability_mass,D_means,D_stds,D_moments)        
        if not self._IsTauleaping:
            self.data_stochsim.setSpecies(sim_dat[:,1:-1].astype(np.uint32),self.sim_species_tracked) 
            self.data_stochsim.setFiredReactions(sim_dat[:,-1].astype(np.uint16))
        else:
            self.data_stochsim.setSpecies(sim_dat[:,1:].astype(np.uint32),self.sim_species_tracked)
            
        self.data_stochsim.setSimulationInfo(self.SSA.timestep,self.SSA.sim_t,self._current_trajectory)
        if self._IsTrackPropensities:
            (L_probability_mass,D_means,D_stds,D_moments) = Analysis.GetDataDistributions(self.SSA.propensities_output,self.sim_rates_tracked)
            self.data_stochsim.setPropensitiesDistributions(L_probability_mass,D_means,D_stds,D_moments)
            self.data_stochsim.setPropensities(self.SSA.propensities_output,self.sim_rates_tracked)
