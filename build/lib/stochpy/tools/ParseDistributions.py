"""
Parse Distributions
===================

Tools to parse model input for delayed, cell division and next reaction methods.

Written by Maxim Moinat and T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: May 26, 2015
"""

import sys, numpy as np

def retrieve_index(input_, lst_names):
    """
    Converts identifiers (or indices) to indices
    
    Input:
     - *name_or_index*
     - *lst_names*  
    """
    if input_ in lst_names:
        index = lst_names.index(input_) # Index according to order in lst_names
    else:
        if isinstance(input_,int) and input_ < len(lst_names): # Treat input as index if this is true
            index = input_
        else:
            raise Warning("'{0}' is not recognized. Choose one of {1} or {2}.".format(input_, lst_names, range(len(lst_names))))            
    return index    
    
def convertInput2Indices(input_,lst_names):
    """
    convert input (indices, identifiers, or a combination) to indices
    
    Input: 
     - *input_* (str,int,list) Example: 'R1',1,['R1',1']
     - *lst_names* (lst) Example: ['R1','R2','R3']
    
    Returns:
       [0], [1], [0,1]
    """    
    output = []
    if type(input_) in [str, int]:
        input_ = [input_]        
    for name in input_:
        index = retrieve_index(name, lst_names)
        output.append(index)    
    return output  

def convertInput2IndicesAndValues(input_, lst_names, default_value):  
    """
    Converts identifiers (or indices) to indices
    
    Input:
     - *input_* (list, dict, tuple)
       Examples: ['R1','R2'] --> [(0,default_value),(1,default_value)]
                 {'R1':2,'R2':4} --> [(0,2),(1,4)]
                 ['R1',('R2',4)] --> [(0,default_value),(1,4)]
                 ['R1',('R1',3)('R2',4)] --> [(0,default_value),(1,4)] Note: we use only the first instance of a reaction
     - *lst_names* (list) of which we extract the indices
     - *default_value* (float) value assigned if no value given
     
    Output:
      - list of (<index>,value) pairs or list of indices, where the index is the place of the name in lst_names
      - or list of indices linking     
    """
    output = []        
    detected_indices = []
    if type(input_) in [str, int]:           # Convert to list if single name or index is specified.
        input_ = [input_]   
              
    for name_value in input_:
        if type(name_value) in [list,set,tuple] and len(name_value) == 2: #The name and value are specified
            index = retrieve_index(name_value[0], lst_names)
            value = name_value[1]
            
        elif type(name_value) in [str, int]: # No value specified, 
            index = retrieve_index(name_value, lst_names)
            if isinstance(input_,dict):
                value = input_[name_value]   # If input_ dictionary, <name_value> is the key.
            else:
                value = default_value                    
        else:
            raise Warning("Unexpected format {0}".format(name_value))
            
        if index not in detected_indices:    # Make sure that each index goes in there only once, i.e. the first one only                  
            output.append((index,value))   
            detected_indices.append(index)   
    return output
    
def ParseDistributions(distributions, reaction_names):
    """
    Returns two dictionaries specifying for each reaction the function and parameter.
    
    Input:
      - *distributions* (dict): dict containing tuples or lists specifying (name, parameters) of the distribution. Example: {"R1": ("fixed",1),'R2':('gamma',0.059,223)}
      - *reaction_names* (list)
     
    Output:
      - *distr_functions*, *distr_parameters* [dictionaries]
      - Keys: index corresponding to list in <reaction_names> 
      - Values: a NumPy random function and list of all parameters, respectively.
      
    #Note: the NumPy exponential distribution uses the SCALE parameter! scale = 1/rate
    """    
    if distributions == None:
        return {},{} #For SMM: no distributions given, then initialize as empty dicts            
    elif type(distributions) == dict:
        distr_functions, distr_parameters = ParseDistributionDictionaries(distributions, reaction_names)
    else:
        raise Warning("Distributions input must be provided as a dictionary")        
    return distr_functions, distr_parameters

def ParseDistributionDictionaries(distr_dict, reaction_names):  
    """ 
    Assigns NumPy distribution functions and parameters to specific reactions specified as keys in input *distr_dict*
    
    We support both reaction identifiers and indices
    
    Input
     - *distr_dict* (dict) Example: {"R1": ("fixed",1),'R2':('gamma',0.059,223)}
     - *reaction_names* (list)    
    """
    functions = {}
    parameters = {}
    for r in distr_dict:        
        if type(r) == str and r in reaction_names:       # r is a string
            r_index = reaction_names.index(r)
        elif type(r) == int and r < len(reaction_names): # r is an index
            r_index = r
        else:
            raise Warning("Reaction name ('{0}') is not found in the model. Choose one of ({1}) or an index.".format(r, reaction_names) )
        
        function, parm = MakeDistributionFunction(distr_dict[r])
        functions[r_index]  = function
        parameters[r_index] = parm
        
    return functions, parameters
    
def MakeDistributionFunction(function_parm):
    """
    Make the function and parameters 
    
    Input: 
     - *function_parm* (tuple or list) distribution name followed by all parameters. Example: ('gamma',0.059,223)
     
    Output:
     - *distr_function* Example:  <function gamma>
     - *distr_parm* Example: (0.059,223)
     
    The distribution NumPy function with the parameters distr_function(*distr_parm, size = 5) generates 5 random variables.    
    """
    distr_name = function_parm[0].lower()
    distr_parm = function_parm[1:]    
    if distr_name == 'fixed':
        distr_function = lambda x, size=None: np.ones(size)*float(x) if size != None else float(x) # allow multiple fixed values produced
        distr_function.__module__ = None  # Purely cosmetic
        distr_function.__name__ = 'fixed' # Purely cosmetic
    else:
        try:
            distr_function = np.random.__dict__[ distr_name ]
        except KeyError:
            raise Warning("The distribution function name '{0}' is not recognized.".format(distr_name) )
    
    #Check whether the distribution function and parameters correspond.
    try:
        float(distr_function(*distr_parm)) #The float() gives error when a list is formed (i.e. the size parameter is assigned)
    except:
        raise Warning("Error: the parameters specified {0} do not correspond to a(n) '{1}' distribution".format(distr_parm, distr_name))     
    return distr_function, distr_parm  

def hypoexp_pdf(x, rates_list):   
    n_rates = len(rates_list)
    
    alpha = np.array([1] + [0]*(n_rates-1))
    theta = np.zeros((n_rates,n_rates))
    for i in range(n_rates-1):
        rate = rates_list[i]
        theta[i][i] = -rate
        theta[i][i+1] = rate
    theta[-1][-1] = -rates_list[-1]
    
    exp_xtheta = linalg.expm(x*theta)
    
    column_vector = np.ones((n_rates,1))
    sumtheta = theta.dot(column_vector)
    
    result = - alpha.dot(exp_xtheta.dot(sumtheta))
    return result[0]
