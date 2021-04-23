 #! /usr/bin/env python
"""
Cain/StochSkit Output to StochPy
================================

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: January 06, 2015
"""
from __future__ import division, print_function, absolute_import

import numpy as np

class Species():
    def __init__(self):
        """ Object that is created to store the species quantities """
        pass
        
def GetCainTimeSim(species_quantities,sim_time,n_frames,n_species): 
    """
    get Cain time simulation output
    
    Input:
     - *species_quantities* (list)
     - *sim_time* (float)
     - *n_frames* (int)
     - *n_species* (int)
    """       
    sim_output = []
    n=0
    species_quantities = [int(s_amount) for s_amount in species_quantities]
    for frame in range(n_frames):            
        time_event = [sim_time[frame]]
        time_event += [species_quantities[n+m] for m in range(n_species)]
        n+=n_species         
        sim_output.append(time_event)   
    return np.array(sim_output,dtype=int)                                        

def GetStochKitTimeSim(file_in,sim_time,species_order):
    """
    get stochkit time simulation output
    
    Input:
     - *file_in* (file)
     - *sim_time* (float)     
     - *species_order* (list)
    """   
    sim_output = []
    IsInit = True
    frame_id = 0   
    for line in file_in:
        dat = line.strip().split('\t')        
        if IsInit and dat[0] == 'time':
            time_event_len = len(dat)
            species = dat[1:]                 
            IsInit = False                
        else:         
            time_event = [int(s_amount) for s_amount in dat[1:]]
            time_event.insert(0,sim_time[frame_id])        
            sim_output.append(time_event)
            frame_id +=1
    #print(species,species_order)    
    Arr_sim_output = np.array(sim_output,dtype=int)  
    if species != species_order:
       Arr_sim_output_cc = np.array(sim_output,dtype=int)  
       for i,s_id in enumerate(species_order):
           s_index = species.index(s_id)
           Arr_sim_output[:,i+1] = Arr_sim_output_cc[:,s_index+1]    
            
    return Arr_sim_output,species_order
        
def GetPropensities(SSAmodule,sim_output):
    """
    get Propensities
    
    Input:
     - *SSAmodule* (python object)
     - *sim_output* (list)
    """
    code_str = """"""
    for i in range(SSAmodule.n_reactions):                    
        code_str += "prop_vec[{0:d}]={1}\n".format(i,SSAmodule.parse.propensities[i])
    req_eval_code = compile(code_str,"RateEqEvaluationCode","exec")
    __species__ = Species()                          
    propensities_output = []
    propensities_distributions = [{} for i in range(SSAmodule.n_reactions)]
    for i in range(len(sim_output)):
        prop_vec = np.zeros(SSAmodule.n_reactions)
        [setattr(__species__,SSAmodule.parse.species[s],sim_output[i][s+1]) for s in range(SSAmodule.n_species)]
        try: 
            exec(req_eval_code)    
        except Exception as er:
            print(er)
        prop_vec = list(prop_vec)            
        prop_vec.insert(0,sim_output[i][0])                  
        propensities_output.append(prop_vec)  
    return propensities_output#,propensities_distributions
