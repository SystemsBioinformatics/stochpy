#! /usr/bin/env python
"""
PySCeS MDL Parser
=================
The PySCeS parser is used to import a model written in the MDL of PySCeS. Further, all required input do to stochastic simulations is build.

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: August 06, 2015
"""

from __future__ import division, print_function, absolute_import

import os,copy,time,sys,operator,re
try: 
    import numpy as np
except ImportError:
    print("Error: The NumPy module is not installed")
    sys.exit()

from stochpy import model_dir, output_dir

from . import PyscesParse
from . import SBML2PSC
from ..core2.InfixParser import MyInfixParser

InfixParser = MyInfixParser()
InfixParser.buildlexer()
InfixParser.buildparser(debug=0, debugfile='infix.dbg', tabmodule='infix_tabmodule',outputdir = output_dir) # 28/08/2014. outputdir added
InfixParser.setNameStr('self.', '')

mach_spec = np.MachAr()
pscParser = PyscesParse.PySCeSParser(debug=0)

class NewCoreBase(object):
    """
    Core2 base class, needed here as we use Core2 derived classes
    in PySCes
    """
    name = None
    __DEBUG__ = False

    def getName(self):
        return self.name

    def setName(self,name):
        self.name = name

    def get(self, attr):
        """Return an attribute whose name is str(attr)"""
        return self.__getattribute__(attr)
        
# this must stay in sync with core2
class NumberBase(NewCoreBase):
    """
    Derived Core2 number class.
    """
    value = None
    value_initial = None

    def __call__(self):
        return self.value

    def getValue(self):
        return self.value

    def setValue(self, v):
        self.value = v
        
class EventAssignment(NumberBase):
    """
    Event assignments are actions that are triggered by an event.
    Ported from Core2 to build an event handling framework fro PySCeS
    """
    variable = None
    symbols = None
    formula = None
    code_string = None
    xcode = None
    mod = None
    piecewises = None
    __DEBUG__ = False

    def __call__(self):
        setattr(self.mod, self.variable, self.value)
        if self.__DEBUG__: print('\tAssigning {0} = {1}'.format(self.variable, self.value))
        return True

    def __init__(self, name, mod):
        self.setName(name)
        self.mod = mod

    def setVariable(self, var):
        self.variable = var

    def setFormula(self, formula):
        self.formula = formula
        InfixParser.setNameStr('self.mod.', '')
        InfixParser.SymbolReplacements = {'_TIME_':'self.mod._TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.code_string = 'self.value={0}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, 'EvAs: {0:s}'.format(self.name), 'exec')
        ##  print( '\t', self.name, self.code_string)

    def evaluateAssignment(self):
        exec(self.xcode)

# adapted from core2
class Event(NewCoreBase):
    """
    Event's have triggers and fire EventAssignments when required.
    Ported from Core2.
    """
    trigger = None
    delay = 0.0
    formula = None
    code_string = None
    xcode = None
    state0 = False
    state = False
    assignments = None
    _TIME_ = 0.0
    _ASS_TIME_ = 0.0   
    _TRIGGER_TIME_ = False
    _need_action = False
    symbols = None
    _time_symbol = None
    piecewises = None
    mod = None
    __DEBUG__ = True

    def __init__(self, name, mod):
        self.setName(name)
        self.assignments = []
        self.mod = mod

    def __call__(self, time,X_matrix):#,__species__):        
        self._TIME_ = time
        exec(self.xcode)               
        ret = False
        if self.state0 and not self.state:
            self.state0 = self.state
        if not self.state0 and self.state:
            #for ass in self.assignments:
            #    ass.evaluateAssignment()
            self.state0 = self.state
            self._need_action = True
            self._ASS_TIME_ = time + self.delay            
            ret = False
        if self._need_action and self._TIME_ >= self._ASS_TIME_:
            #for ass in self.assignments:
            #    ass()
            self._need_action = False
            ret = True
        return ret

    def setTrigger(self, formula, delay=0.0):
        self.formula = formula
        self.delay = delay
        InfixParser.setNameStr('self.mod.', '')
        ##  print formula, delay
        if self._time_symbol != None:
            InfixParser.SymbolReplacements = {self._time_symbol : '_TIME_'}
        else:
            InfixParser.SymbolReplacements = {'_TIME_' : '_TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.code_string = 'self.state={0}'.format(InfixParser.output)
        self.code_string = self.code_string.replace('self.mod._TIME_','self._TIME_')

    def setAssignment(self, var, formula):
        ass = EventAssignment(var, mod=self.mod)
        ass.setVariable(var)
        ass.setFormula(formula)
        self.assignments.append(ass)
        self.__setattr__('_'+var, ass)

    def reset(self):
        self.state0 = False
        self.state = False
        self._TIME_ = 0.0
        self._ASS_TIME_ = 0.0        
        self._TRIGGER_TIME_ = False

class NewCoreBase(object):
    """ Core2 base class, needed here as we use Core2 derived classes in PySCeS """
    name = None
    __DEBUG__ = False

    def getName(self):
        return self.name

    def setName(self,name):
        self.name = name

    def get(self, attr):
        """Return an attribute whose name is str(attr)"""
        return self.__getattribute__(attr)

class Function(NewCoreBase):
    """ Function class ported from Core2 to enable the use of functions in PySCeS """
    formula = None
    code_string = None
    xcode = None
    value = None
    symbols = None
    argsl = None
    mod = None
    piecewises = None
    functions = None

    def __init__(self, name, mod):
        self.name = name
        self.argsl = []
        self.functions = []
        self.mod = mod

    def __call__(self, *args):
        for ar in range(len(args)):
            self.__setattr__(self.argsl[ar], args[ar])
        exec(self.xcode)
        return self.value

    def setArg(self, var, value=None):
        self.__setattr__(var, value)
        self.argsl.append(var)

    def addFormula(self, formula):
        self.formula = formula
        InfixParser.setNameStr('self.', '')
        InfixParser.SymbolReplacements = {'_TIME_':'mod._TIME_'}
        InfixParser.parse(formula)
        self.piecewises = InfixParser.piecewises
        self.symbols = InfixParser.names
        self.functions = InfixParser.functions
        self.code_string = 'self.value={0}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, 'Func: {0:s}'.format(self.name), 'exec')

class PyscesInputFileParser(object):
    """ This class contains the PySCeS model loading """
    ModelDir = None
    ModelFile = None
    ModelOutput = None
    __settings__ = None
    N = None        
    def __init__(self, File, dir, output_dir=None,quiet=False):        
        self.ModelDir = dir
        self.ModelFile = File
        self._IsConverted = False
        if output_dir == None:
            self.ModelOutput = os.getcwd()
        else:
            assert os.path.exists(output_dir),"{0:s} is not a valid path".format(output_dir)
        self.__settings__ = {}
        self.InitialiseInputFile(quiet)    
        self.InitialiseEvents(quiet)           

    def InitialiseInputFile(self,quiet=False):
        """ Parse the input file associated with the PySCeS model instance and assign the basic model attributes """
        self.__parseOK = 1 # check that model has parsed ok
        path_ = os.path.join(self.ModelDir,self.ModelFile)
        try:            
            if os.path.exists(os.path.join(self.ModelDir,self.ModelFile)):
                pass
            else: print('\nInvalid self.ModelFile: ' + os.path.join(self.ModelDir,self.ModelFile) )
        except:
            print('*** WARNING ***: Problem verifying: ' + os.path.join(self.ModelDir,self.ModelFile))
        
        if self.ModelFile[-4:] == '.psc':
            pass
        elif self.ModelFile[-4:].lower() == '.xml':
            #try:
            if not quiet:
                print("Info: extension is .xml")            
            SBML2PSC.SBML2PSC(self.ModelFile,self.ModelDir,quiet=quiet)
            if not quiet:
                print("Info: SBML data is converted into psc data and is stored at: {0:s}".format(model_dir))
            self.ModelFile += '.psc'
            self.ModelDir = model_dir            
            self._IsConverted = True
            #except:
            #    print("Error: Make sure that the libsbml and libxml2 are installed and that the input file is written in SBML format")
            #    print("Info: Use the psc format if both libraries are not available.")
            #    sys.exit()          
        else:
            print('Assuming extension is .psc')
            self.ModelFile += '.psc'

        if not quiet: 
            print('Parsing file: {0:s}'.format(os.path.join(self.ModelDir, self.ModelFile)) )
        pscParser.ParsePSC(self.ModelFile,self.ModelDir,self.ModelOutput,quiet)
        badlist = pscParser.KeywordCheck(pscParser.ReactionIDs)
        badlist = pscParser.KeywordCheck(pscParser.Inits,badlist)

        if len(badlist) != 0:
            print('\n******************************\nPSC input file contains PySCeS keywords please rename them and reload:')
            for item in badlist:
                print('   --> {0}'.format(item) )
            print('******************************\n')
            self.__parseOK = 0

        if self.__parseOK:
            # brett 2008
            self.__nDict__ = pscParser.nDict.copy()
            self.__sDict__ = pscParser.sDict.copy()
            self.__pDict__ = pscParser.pDict.copy()
            self.__uDict__ = pscParser.uDict.copy()            
            self.__eDict__ = pscParser.Events.copy()
            self.__aDict__ = pscParser.AssignmentRules.copy()
            # model attributes are now initialised here brett2008
            self.__InitDict__ = {}
            # set parameters and add to __InitDict__
            for p in list(self.__pDict__):
                setattr(self, self.__pDict__[p]['name'], self.__pDict__[p]['initial'])
                self.__InitDict__.update({self.__pDict__[p]['name'] : self.__pDict__[p]['initial']})
            # set species and add to __InitDict__ and set mod.Xi_init
            for s in list(self.__sDict__):
                setattr(self, self.__sDict__[s]['name'], self.__sDict__[s]['initial'])
                if not self.__sDict__[s]['fixed']:
                    setattr(self, self.__sDict__[s]['name']+'_init', self.__sDict__[s]['initial'])
                self.__InitDict__.update({self.__sDict__[s]['name'] : self.__sDict__[s]['initial']})

            # setup keywords
            self.__KeyWords__ = pscParser.KeyWords.copy()
            if self.__KeyWords__['Modelname'] == None:

                self.__KeyWords__['Modelname'] = self.ModelFile.replace('.psc','')
            if self.__KeyWords__['Description'] == None:
                self.__KeyWords__['Description'] = self.ModelFile.replace('.psc','')

            # if SpeciesTypes undefined assume []
            if self.__KeyWords__['Species_In_Conc'] == None:
                self.__KeyWords__['Species_In_Conc'] = True
            # if OutputType is undefined assume it is the same as SpeciesType
            if self.__KeyWords__['Output_In_Conc'] == None:
                if self.__KeyWords__['Species_In_Conc']:
                    self.__KeyWords__['Output_In_Conc'] = True
                else:
                    self.__KeyWords__['Output_In_Conc'] = False

            # set the species type in sDict according to 'Species_In_Conc'
            for s in list(self.__sDict__):
                if not self.__KeyWords__['Species_In_Conc']:
                    self.__sDict__[s]['isamount'] = True
                else:
                    self.__sDict__[s]['isamount'] = False

            # setup compartments
            self.__compartments__ = pscParser.compartments.copy()
            if len(list(self.__compartments__)) > 0:
                self.__HAS_COMPARTMENTS__ = True
            else:
                self.__HAS_COMPARTMENTS__ = False                
                
            if  self.__aDict__ != {}:
                self.__HAS_ASSIGNMENTS__ = True
                print('Assignment(s) detected.')
            else:  
                self.__HAS_ASSIGNMENTS__ = False

            # no (self.)
            self.__fixed_species__ = copy.copy(pscParser.fixed_species)
            #print "Fixes Species",self.__fixed_species__
            self.__species__ = copy.copy(pscParser.species)
            #print "Species Vector",self.__species__
            self.__parameters__ = copy.copy(pscParser.parameters)
            #print "Parms",self.__parameters__
            self.__reactions__ = copy.copy(pscParser.reactions)
	          #print "Reactions--> rate_eqs",self.__reactions__
            self.__modifiers__ = copy.copy(pscParser.modifiers)
            #print self.__modifiers__
            self.__functions__ = pscParser.Functions.copy() 
            #print self.__functions__
            if self.__functions__ != {}:
                print("*** WARNING ***: StochPy does not support function input")
        else:
            print('\nERROR: model parsing error, please check input file.\n')
        # added in a check for model correctness and human error reporting (1=ok, 0=error)
        if len(pscParser.SymbolErrors) != 0:
            print('\nUndefined symbols:\n{0}'.format(self.SymbolErrors) )
        if not pscParser.ParseOK:
            print('\n\n*****\nModel parsing errors detected in input file '+ self.ModelFile +'\n*****')
            print('\nInput file errors')
            for error in pscParser.LexErrors:
                print(error[0] + 'in line:\t' + str(error[1]) + ' ('+ error[2][:20] +' ...)')
            print('\nParser errors')
            for error in pscParser.ParseErrors:
                try:
                    print(error[0] + '- ' + error[2][:20])
                except:
                    print(error)
            assert pscParser.ParseOK == 1, 'Input File Error'            
        
    def InitialiseEvents(self,quiet=False):
        """ Initialise Events """
        self.__events__ = []
        # for each event
        for e in self.__eDict__:
            ev = Event(e, self)
            ev._time_symbol = self.__eDict__[e]['tsymb']
            ev.setTrigger(self.__eDict__[e]['trigger'], self.__eDict__[e]['delay'])
            # for each assignment
            for ass in self.__eDict__[e]['assignments']:
                ev.setAssignment(ass, self.__eDict__[e]['assignments'][ass])
            self.__events__.append(ev)
            setattr(self, ev.name, ev)
            assert (not '_TIME_' in ev.formula or not ev.delay), "Error: The time event, {0:s}, cannot have a delay".format(ev.name)            
        os.chdir(self.ModelOutput)
        if len(self.__events__) > 0:
            self.__HAS_EVENTS__ = True
            if not quiet:
                print('Event(s) detected.')
        else:
            self.__HAS_EVENTS__ = False

class PySCeS_Connector(PyscesInputFileParser):
    def __init__(self,ModelFile,ModelDir,IsTauleaping = False, IsNRM = False, IsDelayed = False, IsSMM = False,IsQuiet=False):
      """ 
      Use the PySCeS parser to import a set of reactions, which will be used to perform stochastic simulations. Further, some initial stuff that is necessary for stochastic simulations is build.
      
      Input:
       - *ModelFile* (string)
       - *ModelDir* (string)
      """
      self._IsConverted = False
      try:
          self.Mod = PyscesInputFileParser(File = ModelFile, dir = ModelDir,quiet=IsQuiet)
      except Exception as er:
          print(er)
          sys.exit()        
                   
      self.BuildN()      
      if IsDelayed or IsSMM:
          self.BuildN_react_prod()                
      if self.Mod._IsConverted:
          self._IsConverted = True
      self.BuildX()
      self.BuildReactions()        
      self.BuildReactionAffects()
      if IsNRM:        
          self.BuildDependencyGraph()          
      if IsSMM:
          self.BuildSpeciesDependsOn()          
          self.GetReactionOrders()    

      if IsTauleaping:
          self.GetReactionOrders()  
          self.species_HORs = np.zeros(len(self.species)).tolist()      
          for i,s_id in enumerate(self.species):
              for j,r_reactants in enumerate(self.reactants):
                  if s_id in r_reactants:
                      if self.reaction_orders[j] > self.species_HORs[i]:
                          self.species_HORs[i] = self.reaction_orders[j]          

    def GetReactionOrders(self):
        """
        Get reaction orders for each reaction. Note that we assume mass-action kinetics
               
        1. Use "DependsOn"  
      
        """
        self.reaction_orders = []        
        self.species_max_influence = np.zeros(len(self.species)).tolist() # species with highest influence [('S1',2),('S2',1)]
        self.reaction_orders = [len(x) for x in self.depends_on]
        for i in range(len(self.species)):
            s_max_influence = max([x.count(i) for x in self.depends_on])
            self.species_max_influence[i] = s_max_influence              

            
    def BuildN(self):
        """ Generates the stoichiometric matrix N from the parsed model description. Returns a stoichiometric matrix (N) as a NumPy array """
        species_self = ['self.'+s for s in self.Mod.__species__] # species identifiers with "self."
        self.N_matrix = np.zeros((len(species_self),len(self.Mod.__reactions__))) #,dtype=np.int16)
        for s_id in species_self:
            for r_id in self.Mod.__reactions__:
                r_reagents = self.Mod.__nDict__[r_id]['Reagents']
                if s_id in r_reagents:
                    s_index = species_self.index(s_id)
                    stoichiometry = r_reagents[s_id]
                    self.N_matrix[s_index][self.Mod.__reactions__.index(r_id)] = stoichiometry

    def BuildN_react_prod(self): #24-10-2013, used in delayed methods
        """
        Generates the reactant, product and net stoichiometric matrix N from the parsed model description. 
        Also: generates list of reactant and product indices for each reaction.
        Returns three stoichiometric matrices (N) as NumPy arrays
        """
        species_self = ['self.'+s for s in self.Mod.__species__] # used for the correct order of species
        self.N_matrix_reactants = np.zeros((len(species_self),len(self.Mod.__reactions__)))
        self.N_matrix_products = np.zeros((len(species_self),len(self.Mod.__reactions__)))
        
        self.reactant_indices = [] # Per reaction the reactant indices 
        self.product_indices = []  # Per reaction the product indices         
        for r_index, r_id in enumerate(self.Mod.__reactions__):
            all_reagents_stoich = self.Mod.__nDict__[r_id]['AllReagents'] #list of tuples with species_name and stoichiometry (N)
            temp_reactants = [] #29-10-2013 -> For SMM we want double occurances (i.e. two products formed). WARNING: This leads to double iterations in delayed method (updating species).
            temp_products = []
            for s_name, s_stoich in all_reagents_stoich:             
                if s_name in species_self:
                    s_index = species_self.index(s_name)
                    s_stoich = int(s_stoich)
                    if s_stoich<0:   # reactant
                        self.N_matrix_reactants[s_index][r_index] += s_stoich
                        temp_reactants.extend( [s_index]*abs(s_stoich) ) #Make sure that length corresponds with number of molecules formed (for SMM)
                    elif s_stoich>0: # product
                        self.N_matrix_products[s_index][r_index] += s_stoich
                        temp_products.extend( [s_index]*s_stoich ) 
            self.reactant_indices.append( list(temp_reactants) )
            self.product_indices.append(  list(temp_products) )
   
    def BuildX(self):  
      """ Builds the initial concentrations of all species (X). """
      self.species = copy.deepcopy(self.Mod.__species__)	   # Species names from parser
      n_species = len(self.species)      
      self.fixed_species_amount = []
      for species in self.Mod.__fixed_species__:
          s_amount = self.Mod.__sDict__[species]['initial']          
          self.fixed_species_amount.append(s_amount)
     
      self.X_matrix = np.zeros((n_species,1)) #,dtype = np.int16) 
      for i,s_id in enumerate(self.species):          
          init_amount = self.Mod.__sDict__[s_id]['initial']
          if not self.Mod.__sDict__[s_id]['isamount']:    # Handle amount/concentration 19 September 2011
              if self.Mod.__compartments__:
                  L_compartments = list(self.Mod.__compartments__)
                  if len(L_compartments) > 1:
                      print("Info: Multiple compartments are detected")
                  for compartment in L_compartments:
                      if self.Mod.__compartments__[compartment]['size'] != 1:
                          if self.Mod.__sDict__[s_id]['compartment'] == compartment:
                              print("*** WARNING ***: Species are given in concentrations and the compartment volume is unequal to one.\nThe species concentrations are multiplied by the compartment volume.")
                              init_amount *= self.Mod.__compartments__[compartment]['size']
          self.X_matrix[i] = init_amount                       # Add initial amount

    def BuildReactions(self):
        """
        Extract information out of each reaction, such as what are the reagents/reactants and which parameter is used for that particular reaction. 
        
        Reactants(r): set of reactants of reaction r
        Products(r): set of products of reaction r
        DependsOn(propensity): set of substances that affect the propensity
        Affects(r): set of substances that change quantity when reaction r is executed        
        """        
        self.reactants = [] 
        self.depends_on = []          
        #self.rate_eqs  = []
        self.propensities = []   
        regex = "self.(\w+)"  
        
        IsReversible = False
        for r_id in self.Mod.__reactions__:                    
            rate_eq = self.Mod.__nDict__[r_id]['RateEq']     # Rate eq. info                 
            temp_depends_on = re.findall(regex,rate_eq)      # parameters and reactants/reagents
            
            for parm in self.Mod.__nDict__[r_id]['Params']:  # Replace parm values (fixed)                  
                parm = parm.replace('self.','')  
                while parm in temp_depends_on:               # remove parameters from depends_on
                    temp_depends_on.remove(parm)                
                if parm in self.Mod.__parameters__:
                    if parm in self.Mod.__pDict__:
                        parm_value = self.Mod.__pDict__[parm]['initial']             
                    elif parm in self.Mod.__sDict__:         # Fix 28 Feb 2013, however, not sure if this is still possible. Let's keep it here just to be sure
                        parm_value = self.Mod.__sDict__[parm]['initial']
                elif parm in self.Mod.__compartments__:
                    parm_value = self.Mod.__compartments__[parm]['size']
                if parm not in self.Mod.__fixed_species__:   # March 26, ignore fixed species
                    rate_eq = rate_eq.replace('self.{0:s}'.format(parm), str(parm_value))            

            self.depends_on.append([self.species.index(s_id) for s_id in temp_depends_on])
                        
            rate_eq = rate_eq.replace('self.','__species__.') # Make a special object 'species' where all species amounts are stored 20 Oct 2011
            self.propensities.append(rate_eq)                 # Add rate equation to props eq         
            #self.rate_eqs.append(self.Mod.__nDict__[r_id])   #self.rate_eqs.append(self.Mod.__nDict__[r_id]['RateEq'])          

            temp_reactants = []
            for reagent in self.Mod.__nDict__[r_id]['AllReagents']:
                reagent_id = reagent[0].replace('self.','')
                if reagent[-1] < 0 and reagent_id not in temp_reactants:                
                    temp_reactants.append(reagent_id)
             
            self.reactants.append(temp_reactants)      
            if self.Mod.__nDict__[r_id]['Type'] == 'Rever':        
                IsReversible = True
            assert not IsReversible, "Error: The model contains reversible reactions, while Stochastic Simulation Algorithms require irreversible reactions."

    def BuildReactionAffects(self):  
        """
        Determine the affects for each reaction
        
        Affects(r): set of substances that change quantity when reaction r is executed                
        """
        self.reaction_affects = []          
        for r_id in self.Mod.__reactions__: 
            D_reagents = self.Mod.__nDict__[r_id]['Reagents']            
            affects_r = []      
            for species in D_reagents:               
                s_id = species.replace('self.','')  # September 19, 2013    
                if s_id in self.species: # only variable species
                    affects_r.append(self.species.index(s_id))            
            self.reaction_affects.append(affects_r)

    def BuildDependencyGraph(self):
        """ Function which builds a dependency graph """
        # Bug removed on November 06, 2013       
        self.dep_graph = []       
        for i in range(len(self.reaction_affects)):
            to_change = []
            for s_index in self.reaction_affects[i]:                  
                for j,r_indices in enumerate(self.depends_on):
                    if s_index in r_indices:
                       to_change.append(j)                    
            to_change.append(i) # 21-11-2013: add own reaction index.
            self.dep_graph.append(np.unique(to_change).tolist())
      
      
    def BuildSpeciesDependsOn(self): # 29-10-2013 -> used in Single Molecule Method
        """        
        Builds dependencies specifying for each species in which reactions they are reactants and/or modifiers.
        Used for the SMM method(s).
        """
        self.species_depends_on = []
        for s_index, s_id in enumerate(self.species):    
            change = set()                        # Do not allow duplicates of a reaction             
            for j, reactants in enumerate(self.depends_on): # 28-04-2014, replaces enumerate(self.reactant_indices): 
                if s_index in reactants:          # Species is reactant of reaction j
                    change.add(j)
                
                elif s_id in self.reactants[j]:   # Backup for modifiers, if they are not a reactant.
                    print("*** WARNING ***: rate equation of reaction '{0:s}' contains a species which is not a reactant.".format(self.Mod.__reactions__[j]) )
                    change.add(j)                    
            self.species_depends_on.append( list(change) )    
      
      
class RegularGridDataObj(object):    
    """
    This class is specifically designed to store the results of a stochastic time simulation on fixed time intervals
    It has methods for setting the Time, Labels, Species and Propensity data and
    getting Time, Species and Rate (including time) arrays. However, of more use:

    - getOutput(\*args) feed this method species/rate labels and it will return
      an array of [time, sp1, r1, ....]
    - getDataAtTime(time) the data generated at time point "time".
    - getDataInTimeInterval(time, bounds=None) more intelligent version of the above
      returns an array of all data points where: time-bounds <= time <= time+bounds
    """
    time_label = 'Time'    
    time = None
    species = None
    propensities = None
    species_autocorrelations = None 

    species_means = None
    species_standard_deviations = None
    species_autocorrelations_means = None
    species_autocorrelations_standard_deviations = None
    species_labels = None
    propensities_means = None
    propensities_standard_deviations = None    
    propensities_labels = None    
    
    HAS_SPECIES = False
    HAS_VOLUME = False
    HAS_SPECIES_CONCENTRATIONS = False
    HAS_PROPENSITIES = False    
    HAS_TIME = False  
    HAS_SPECIES_AUTOCORRELATIONS = False
    HAS_PROPENSITIES_AUTOCORRELATIONS = False
    
    HAS_SPECIES_AUTOCOVARIANCES = False
    HAS_PROPENSITIES_AUTOCOVARIANCES = False
    
    HAS_AVERAGE_SPECIES_DISTRIBUTIONS = False
    HAS_AVERAGE_PROPENSITIES_DISTRIBUTIONS = False    
    
    def setSpeciesDistributionAverage(self,mean,std):
        """
        Set means and stds of species data
        
        Input:
         - *mean* (list)
         - *std* (list)
        """
        self.species_distributions_means = mean
        self.species_distributions_standard_deviations = std
        self.HAS_AVERAGE_SPECIES_DISTRIBUTIONS = True

    def setPropensitiesDistributionAverage(self,mean,std):
        """
        Set means and stds of species data

        Input:
         - *mean* (list)
         - *std* (list)
        """
        self.propensities_distributions_means = mean
        self.propensities_distributions_standard_deviations = std        
        self.HAS_AVERAGE_PROPENSITIES_DISTRIBUTIONS = True
    
    def setSpecies(self, species, lbls=None):
        """
        Set the species array
        
        Input:
         - *species* an array of species vs time data
         - *lbls* [default=None] a list of species labels
        """
        self.species = species
        self.HAS_SPECIES = True
        if lbls != None:
            self.species_labels = lbls   
            
    def setVolume(self,volume):
        """
        Set the volume array
        
        Input:
         - *species* an array of volume vs time data
        """
        self.volume = np.array(volume)
        self.HAS_VOLUME = True          

    def setSpeciesConcentrations(self,extracellular=[]):
        """ Set the species concentrations array """
        assert self.HAS_SPECIES,"Error: species not set yet"
        assert self.HAS_VOLUME, "Error: volume not set yet"

        self.species_concentrations = self.species/self.volume
        self.species_concentrations_labels = self.species_labels        
        if extracellular != []:               
            self.species_concentrations = np.delete(self.species_concentrations,extracellular,1)
            self.species_concentrations_labels = np.delete(self.species_concentrations_labels,extracellular,0).tolist()
            
        self.HAS_SPECIES_CONCENTRATIONS = True        
            
    def setTime(self, time, lbl=None):
        """
        Set the time vector

        Input:
         - *time* a 1d array of time points
         - *lbl* [default=None] is "Time" set as required
        """
        self.time = time
        self.HAS_TIME = True
        if lbl != None:
            self.time_label = lbl

    def getSpecies(self, lbls=False):
        """
        Return an array of time+species

        Input:
        - *lbls* [default=False] return only the time+species array or optionally both the data array and a list of column label
        """
        output = None
        if self.HAS_SPECIES:
            output = np.column_stack((self.time, self.species_means))
            labels = [self.time_label]+self.species_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels  
            
    def getPropensities(self, lbls=False):
        """
        Return time+propensity array

        Input:        
         - *lbls* [default=False] return only the time+propensity array or optionally both the data array and a list of column label        

        """
        assert self.propensities != None, "\nNo propensities"
        output = None
        if self.HAS_PROPENSITIES:
            output = np.column_stack((self.time, self.propensities_means))
            labels = [self.time_label]+self.propensities_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels            
            
            
    def getVolume(self):
        """ Return an array of time+volume """
        output = None
        if self.HAS_VOLUME:
            output = np.column_stack((self.time, self.volume))           
        else:
            output = self.time              
        return output        
          
            
    def getSpeciesConcentrations(self, lbls=False):
        """
        Return an array of time+species concentrations

        Input:
        - *lbls* [default=False] return only the time+species array or optionally both the data array and a list of column label
        """
        output = None
        if self.HAS_SPECIES_CONCENTRATIONS:
            output = np.column_stack((self.time, self.species_concentrations))
            labels = [self.time_label]+self.species_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels               

            
    def getTime(self, lbls=False):
        """
        Return the time vector

        Input:
         - *lbls* [default=False] return only the time array or optionally both the time array and time label
        """
        output = None
        if self.HAS_TIME:
            output = self.time
        if not lbls:
            return output
        else:
            return output, [self.time_label]    
            
    def setPropensities(self, propensities,lbls=None):
        """
        Sets an array of propensities.

        Input:
         - *propensities* (list)
        """ 
        self.propensities = propensities
        self.HAS_PROPENSITIES = True 
        if lbls != None:
            self.propensities_labels = lbls            
        
    def setPropensitiesLabels(self,labels): 
        """
        Input:
         - *labels* (list)
        """
        self.propensities_labels = labels                            
        

    def setSpeciesAutocorrelations(self,auto_correlations, lbls=None):
        """
        Set the `autocorrelations` ***

        Input:
         - *auto_correlations* (list)
         - *lbls* [default=None] a list of matching reaction names         
        """
        self.species_autocorrelations = np.array(auto_correlations)
        self.HAS_SPECIES_AUTOCORRELATIONS = True   

    def setSpeciesAutocovariances(self,auto_covariances, lbls=None):
        """
        Set the `autocorrelations` ***


        Input:
         - *auto_covariances* (list)
         - *lbls* [default=None] a list of matching reaction names         
        """
        self.species_autocovariances = np.array(auto_covariances)
        self.HAS_SPECIES_AUTOCOVARIANCES = True       
        
    def setPropensitiesAutocorrelations(self,auto_correlations, lbls=None):
        """
        Set the `autocorrelations` ***

        Input:
         - *auto_correlations* (list)
         - *lbls* [default=None] a list of matching reaction names         
        """
        self.propensities_autocorrelations = np.array(auto_correlations)
        self.HAS_PROPENSITIES_AUTOCORRELATIONS = True

    def setPropensitiesAutocovariances(self,auto_covariances, lbls=None):
        """
        Set the `autocovariances` ***


        Input:
         - *auto_covariances* (list)
         - *lbls* [default=None] a list of matching reaction names         
        """
        self.propensities_autocovariances = np.array(auto_covariances)
        self.HAS_PROPENSITIES_AUTOCOVARIANCES = True        
        
             

class IntegrationStochasticDataObj(object):
    """
    This class is specifically designed to store the results of a stochastic time simulation
    It has methods for setting the Time, Labels, Species and Propensity data and
    getting Time, Species and Rate (including time) arrays. However, of more use:

    - getOutput(\*args) feed this method species/rate labels and it will return
      an array of [time, sp1, r1, ....]
    - getDataAtTime(time) the data generated at time point "time".
    - getDataInTimeInterval(time, bounds=None) more intelligent version of the above
      returns an array of all data points where: time-bounds <= time <= time+bounds
    """
    time = None
    waiting_times = None
    species = None
    species_distributions = None       
    
    propensities = None
    propensities_distributions = None
    
    xdata = None
    time_label = 'Time'
    waiting_times_labels = None
    species_labels = None
    propensities_labels = None
    xdata_labels = None
    HAS_SPECIES = False
    HAS_SPECIES_CONCENTRATIONS = False
    HAS_VOLUME = False
    HAS_WAITINGTIMES = False
    HAS_PROPENSITIES = False    
    HAS_TIME = False
    HAS_XDATA = False
    IS_VALID = True
    TYPE_INFO = 'Stochastic'
    
    def setSpeciesDistributions(self,distributions,means,stds,moments):
        """
        setSpeciesDistributions() stuff for the determination of distributions
        
        Input:
         - *distributions* (list)
         - *means* (dictionary)
         - *stds* (dictionary)
         - *moments* (dictionary) 
        """
        self.species_distributions = distributions
        self.species_means = means
        self.species_standard_deviations = stds
        self.species_moments = moments
        

    def setPropensitiesDistributions(self,distributions,means,stds,moments):
        """
        setPropensitiesDist stuff for the determination of distributions
        
        Input:
         - *distributions* (list)
         - *means* (dictionary)
         - *stds* (dictionary)
         - *moments* (dictionary) 
        """
        self.propensities_distributions = distributions
        self.propensities_means = means
        self.propensities_standard_deviations = stds
        self.propensities_moments = moments
        

    def setSimulationInfo(self,timesteps,endtime,simulation_trajectory):
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
        
        
    def setFiredReactions(self,fired_reactions):
        """
        Set the reactions that fired
        
        Input:    
         - *fired_reactions* (list)
        """
        self.fired_reactions = fired_reactions    
            

    def setLabels(self, species):
        """
        Set the species
        
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
        self.time = time
        self.HAS_TIME = True
        if lbl != None:
            self.time_label = lbl
            

    def setSpecies(self, species, lbls=None):
        """
        Set the species array
        
        Input:
         - *species* an array of species vs time data
         - *lbls* [default=None] a list of species labels
        """
        self.species = species
        self.HAS_SPECIES = True
        if lbls != None:
            self.species_labels = lbls
            

    def setVolume(self,volume):
        """
        Set the volume array
        
        Input:
         - *species* an array of volume vs time data
        """
        self.volume = np.array(volume)
        self.HAS_VOLUME = True      


    def setVolumeDistribution(self,distribution):
        """
        set volume distribution
        
        Input:
         - *distribution* (list)
        """
        self.volume_distribution = distribution
                 
            
    def setSpeciesConcentrations(self,extracellular=[]):
        """ Set the species concentrations array """
        assert self.HAS_SPECIES,"Error: species not set yet"
        assert self.HAS_VOLUME, "Error: volume not set yet"        
          
        self.species_concentrations = self.species/self.volume.reshape(len(self.volume), 1)
        self.species_concentrations_labels = self.species_labels
        if extracellular != []:               
            self.species_concentrations = np.delete(self.species_concentrations,extracellular,1)
            self.species_concentrations_labels = np.delete(self.species_concentrations_labels,extracellular,0).tolist()
            
        self.HAS_SPECIES_CONCENTRATIONS = True       
        
            
    def getSpeciesConcentrations(self, lbls=False):
        """
        Return an array of time+species concentrations

        Input:
        - *lbls* [default=False] return only the time+species array or optionally both the data array and a list of column label
        """
        output = None
        if self.HAS_SPECIES_CONCENTRATIONS:
            output = np.column_stack((self.time, self.species_concentrations))
            labels = [self.time_label]+self.species_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels              
            

    def setWaitingtimesMeans(self,waiting_times,rate_names):
        """
        set waiting times means
        
        Input:
         - *waiting_times* (dictionary)
        """
        self.waiting_times_means = []
        for r_id in rate_names:
            waiting_times_r = waiting_times[r_id]
            self.waiting_times_means.append(np.mean(waiting_times_r))
            

    def setWaitingtimesStandardDeviations(self,waiting_times,rate_names):
        """
        set waiting times standard deviations
        
        Input:
         - *waiting_times* (dictionary)
        """
        self.waiting_times_standard_deviations = []       
        for i,r_id in enumerate(rate_names):
            waiting_times_r = np.array(waiting_times[r_id])
            if len(waiting_times_r):                
                variance = sum((self.waiting_times_means[i]- waiting_times_r)**2)/len(waiting_times_r)
                self.waiting_times_standard_deviations.append(variance**0.5)         
            else:
                self.waiting_times_standard_deviations.append(np.NAN)            
          

    def setWaitingtimes(self, waiting_times, lbls=None):
        """
        Set the `waiting_times` this data structure is not an array but a nested list of: waiting time log bins per reaction per trajectory.         
        waiting_times = [traj_1, ..., traj_n]
        traj_1 = [wt_J1, ..., wt_Jn] # in order of SSA_REACTIONS
        wt_J1 = (xval, yval, nbin)
        xval =[x_1, ..., x_n]
        yval =[y_1, ..., y_n]
        nbin = n

        Input:
         - *waiting_times* a list of waiting times
         - *lbls* [default=None] a list of matching reaction names         
        """
        self.waiting_times = waiting_times
        self.HAS_WAITINGTIMES = True
        if lbls != None:
            self.waiting_times_labels = lbls
            

    def setPropensitiesLabels(self,labels): 
        """
        Input:
         - *labels* (list)
        """
        self.propensities_labels = labels     
               

    def setPropensities(self, propensities,lbls=None):
        """
        Sets an array of propensities.

        Input:
         - *propensities* (list)
        """ 
        P_ARR = np.zeros((len(propensities), len(propensities[0])-1), 'd')        
        #P_ARR[-1,:] = np.NaN
        for r in range(P_ARR.shape[0]):            
            P_ARR[r, :] = propensities[r][1:]
        self.propensities = P_ARR
        self.HAS_PROPENSITIES = True 
        if lbls != None:
            self.propensities_labels = lbls  
                     

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
            

    def getTime(self, lbls=False):
        """
        Return the time vector

        Input:
         - *lbls* [default=False] return only the time array or optionally both the time array and time label
        """
        output = None
        if self.HAS_TIME:
            output = self.time
        if not lbls:
            return output
        else:
            return output, [self.time_label]
            

    def getSpecies(self, lbls=False):
        """
        Return an array of time+species

        Input:
        - *lbls* [default=False] return only the time+species array or optionally both the data array and a list of column label
        """
        output = None
        if self.HAS_SPECIES:
            output = np.column_stack((self.time, self.species))
            labels = [self.time_label]+self.species_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels
            
            
    def getVolume(self):
        """ Return an array of time+volume """
        output = None
        if self.HAS_VOLUME:
            output = np.column_stack((self.time, self.volume))
        else:
            output = self.time            
        return output             
        
        
    def getPropensities(self, lbls=False):
        """
        Return time+propensity array

        Input:        
         - *lbls* [default=False] return only the time+propensity array or optionally both the data array and a list of column label        
        """
        assert self.propensities != None, "\nNo propensities"
        output = None
        if self.HAS_PROPENSITIES:
            output = np.column_stack((self.time, self.propensities))
            labels = [self.time_label]+self.propensities_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels
            

    def getXData(self, lbls=False):
        """
        Return time+xdata array

        Input:
        - *lbls* [default=False] return only the time+xdata array or optionally both the data array and a list of column label
        """
        output = None
        if self.HAS_XDATA:
            output = np.column_stack((self.time, self.xdata))
            labels = [self.time_label]+self.xdata_labels
        else:
            output = self.time
            labels = [self.time_label]
        if not lbls:
            return output
        else:
            return output, labels
            

    def getDataAtTime(self, time):
        """
        Return all data generated at "time"

        Input:
         - *time* the required exact time point
        """
        #TODO add rate rule data
        t = None
        sp = None
        ra = None
        ru = None
        xd = None
        
        for tt in range(self.simulation_timesteps+1):
            if self.time[tt] == time:
                t = tt
                if self.HAS_SPECIES:
                    sp = self.species.take([tt], axis=0)
                if self.HAS_PROPENSITIES:
                    ru = self.propensities.take([tt], axis=0)
                if self.HAS_XDATA:
                    xd = self.xdata.take([tt], axis=0)
                break

        output = None
        if t is not None:
            output = np.array([[self.time[t]]])
            if sp is not None:
                output = np.hstack((output,sp)) 
            if ra is not None:
                output = np.hstack((output,ra))
            if ru is not None:
                output = np.hstack((output,ru))
            if xd is not None:
                output = np.hstack((output,xd))
        return output


    def getDataInTimeInterval(self, time, bounds=None):
        """
        Returns an array of all data in interval: time-bounds <= time <= time+bounds
        where bound defaults to step size

        Input:
         - *time* the interval midpoint
         - *bounds* [default=None] interval half span defaults to step size
        """       
        if bounds == None:
            bounds = self.time[1] - self.time[0]
        c1 = (self.time >= time-bounds)
        c2 = (self.time <= time+bounds)
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
            if self.HAS_SPECIES and self.HAS_TIME:
                output = np.hstack((output, self.species.take(t, axis=0)))
            if self.HAS_PROPENSITIES:
                output = np.hstack((output, self.propensities.take(t, axis=0)))
            if self.HAS_XDATA:
                output = np.hstack((output, self.xdata.take(t, axis=0)))
        return output
        

    def getAllSimData(self,lbls=False):
        """
        Return an array of time + all available simulation data

        Input:
         - *lbls* [default=False] return only the data array or (data array, list of labels)
        """
        labels = [self.time_label]
        if self.HAS_SPECIES and self.HAS_TIME:
            output = np.column_stack((self.time, self.species))
            labels += self.species_labels
        if self.HAS_PROPENSITIES:
            output = np.column_stack((output, self.propensities))
            labels += self.propensities_labels
        if self.HAS_XDATA:
            output = np.column_stack((output, self.xdata))
            labels += self.xdata_labels
        if not lbls:
            return output
        else:
            return output, labels
            

    def getSimData(self, *args, **kwargs):
        """
        Feed this method species/xdata labels and it will return an array of [time, sp1, ....]
        
        Input:
         - 'speces_l', 'xdatal' ...
         - *lbls* [default=False] return only the data array or (data array, list of labels)
         
        Example:
        getSimData('mRNA','R1','R2') for the immigration-death model
        """
        output = self.time

        if kwargs.has_key('lbls'):
            lbls = kwargs['lbls']
        else:
            lbls = False
        lout = [self.time_label]
        for roc in args:
            if self.HAS_SPECIES and roc in self.species_labels:
                lout.append(roc) 
                output = np.column_stack((output, self.species.take([self.species_labels.index(roc)], axis=-1)))
            if self.HAS_PROPENSITIES and roc in self.propensities_labels:
                lout.append(roc)
                output = np.column_stack((output, self.propensities.take([self.propensities_labels.index(roc)], axis=-1)))
            if self.HAS_XDATA and roc in self.xdata_labels:
                lout.append(roc)
                output = np.column_stack((output, self.xdata.take([self.xdata_labels.index(roc)], axis=-1)))
        if not lbls:
            return output
        else:
            return output, lout
