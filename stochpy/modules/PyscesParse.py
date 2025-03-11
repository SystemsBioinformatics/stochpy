"""

Copyright (C) 2004-2014 B.G. Olivier

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa.

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license. See LICENSE.txt that came with
this distribution for specifics.
pys
NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier

Adapted by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: February 11, 2015
"""
from __future__ import division, print_function, absolute_import

from stochpy import model_dir, output_dir

import os,copy,sys,importlib
from ..lib import lex
from ..lib import yacc

from getpass import getuser
from time import sleep, strftime

RESERVEDTERMS = (\
     'Vvec',\
     '__CDvarUpString__',\
     '__SALL__',\
     '__SI__',\
     '__epmatrix__',\
     '__evmatrix__',\
     '__fixed_species__',\
     '__kmatrix__',\
     '__kzeromatrix__',\
     '__lmatrix__',\
     '__lzeromatrix__',\
     '__mapFunc_R__',\
     '__mapFunc__',\
     '__mode_suppress_info__',\
     '__modifiers__',\
     '__nmatrix__',\
     '__nrmatrix__',\
     '__parameters__',\
     '__reactions__',\
     '__s_varFunc__',\
     '__species__',\
     '__vFunc__',\
     'cc_all',\
     'cc_all_col',\
     'cc_all_row',\
     'cc_conc',\
     'cc_conc_col',\
     'cc_conc_row',\
     'cc_flux',\
     'cc_flux_col',\
     'cc_flux_row',\
     'conservation_matrix',\
     'conservation_matrix_col',\
     'conservation_matrix_row',\
     'display_debug',\
     'eModeExe_dbl',\
     'eModeExe_int',\
     'elas_epar_upsymb',\
     'elas_evar_upsymb',\
     'elas_par',\
     'elas_par_col',\
     'elas_par_row',\
     'elas_par_u',\
     'elas_var',\
     'elas_var_col',\
     'elas_var_row',\
     'elas_var_u',\
     'emode_file',\
     'emode_intmode',\
     'emode_userout',\
     'fintslv_range',\
     'fintslv_rmult',\
     'fintslv_step',\
     'fintslv_tol',\
     'fixed_species',\
     'hybrd_epsfcn',\
     'hybrd_factor',\
     'hybrd_maxfev',\
     'hybrd_mesg',\
     'hybrd_xtol',\
     'info_flux_conserve',\
     'info_moiety_conserve',\
     'init_species_array',\
     'kel_unscaled',\
     'kmatrix',\
     'kmatrix_col',\
     'kmatrix_row',\
     'kmatrix_scaled',\
     'kzeromatrix',\
     'kzeromatrix_col',\
     'kzeromatrix_row',\
     'lmatrix',\
     'lmatrix_col',\
     'lmatrix_row',\
     'lmatrix_scaled',\
     'lsoda_atol',\
     'lsoda_h0',\
     'lsoda_hmax',\
     'lsoda_hmin',\
     'lsoda_mesg',\
     'lsoda_mxordn',\
     'lsoda_mxords',\
     'lsoda_mxstep',\
     'lsoda_rtol',\
     'lzeromatrix',\
     'lzeromatrix_col',\
     'lzeromatrix_row',\
     'mach_floateps',\
     'mca_ccall_altout',\
     'mca_ccall_concout',\
     'mca_ccall_fluxout',\
     'mca_ccj_upsymb',\
     'mca_ccs_upsymb',\
     'mca_ci',\
     'mca_ci_col',\
     'mca_ci_row',\
     'mca_cjd',\
     'mca_cjd_col',\
     'mca_cjd_row',\
     'mca_csd',\
     'mca_csd_col',\
     'mca_csd_row',\
     'misc_write_arr_lflush',\
     'mode_eigen_output',\
     'mode_elas_deriv',\
     'mode_elas_deriv_dx',\
     'mode_elas_deriv_factor',\
     'mode_elas_deriv_min',\
     'mode_elas_deriv_order',\
     'mode_mca_scaled',\
     'mode_number_format',\
     'mode_sim_init',\
     'mode_sim_max_iter',\
     'mode_solver',\
     'mode_solver_fallback',\
     'mode_solver_fallback_integration',\
     'mode_state_init',\
     'mode_state_init2_array',\
     'mode_state_init3_factor',\
     'mode_state_mesg',\
     'modifiers',\
     'nleq2_ibdamp',\
     'nleq2_iormon',\
     'nleq2_iscal',\
     'nleq2_iter',\
     'nleq2_jacgen',\
     'nleq2_mesg',\
     'nleq2_mprerr',\
     'nleq2_nonlin',\
     'nleq2_qnscal',\
     'nleq2_qrank1',\
     'nleq2_rtol',\
     'nmatrix',\
     'nmatrix_col',\
     'nmatrix_row',\
     'nrmatrix',\
     'nrmatrix_col',\
     'nrmatrix_row',\
     'parameters',\
     'pitcon_abs_tol',\
     'pitcon_allow_badstate',\
     'pitcon_filter_neg',\
     'pitcon_filter_neg_res',\
     'pitcon_fix_small',\
     'pitcon_flux_gen',\
     'pitcon_init_par',\
     'pitcon_iter',\
     'pitcon_jac_opt',\
     'pitcon_jac_upd',\
     'pitcon_limit_point_idx',\
     'pitcon_limit_points',\
     'pitcon_max_grow',\
     'pitcon_max_step',\
     'pitcon_max_steps',\
     'pitcon_min_step',\
     'pitcon_output_lvl',\
     'pitcon_par_opt',\
     'pitcon_par_space',\
     'pitcon_rel_tol',\
     'pitcon_start_dir',\
     'pitcon_start_step',\
     'pitcon_targ_val',\
     'pitcon_targ_val_idx',\
     'pitcon_target_points',\
     'reactions',\
     'scan1_dropbad',\
     'scan1_mca_mode',\
     'scan_in',\
     'scan_out',\
     'scan_res',\
     'sim_end',\
     'sim_points',\
     'sim_res',\
     'sim_start',\
     'sim_time',\
     'species',\
     'state_flux',\
     'state_set_conserve',\
     'state_species',\
     'write_array_header',\
     'write_array_html_footer',\
     'write_array_html_format',\
     'write_array_html_header',\
     'write_array_spacer',\
     'zero_val',\
    )

class PySCeSParser:
    """The original PySCeS parser has been rewritten and extended to include
    many features that allow for SBML compatibility such as MathML 2.0 support,
    function definitions, assignment/rate rules and events."""

    MathmlToNumpy_funcs = {
        'pow' : 'pow', 'root' : 'pow', 'abs' : 'abs',
        'exp' : 'math.exp', 'ln' : 'math.log', 'log' : 'math.log10',
        'floor' : 'numpy.floor', 'ceiling' : 'numpy.ceil', 'factorial' : None,
        'sin' : 'numpy.sin', 'cos' : 'numpy.cos', 'tan' : 'numpy.tan',
        'sec' : None, 'csc' : None, 'cot' : None,
        'sinh' : 'numpy.sinh', 'cosh' : 'numpy.cosh','tanh' : 'numpy.tanh',
        'sech' : None, 'csch' : None, 'coth' : None,
        'arcsin' : 'numpy.arcsin', 'arccos' : 'numpy.arccos', 'arctan' : 'numpy.arctan',
        'arcsec' : None, 'arccsc' : None, 'arccot' : None,
        'arcsinh' : 'numpy.arcsinh', 'arccosh' : 'numpy.arccosh', 'arctanh' : 'numpy.arctanh',
        'arcsech' : None, 'arccsch' : None, 'arccoth' : None,
        'eq' : 'operator.eq', 'neq' : 'operator.ne',
        'gt' : 'operator.gt', 'geq' : 'operator.ge',
        'lt' : 'operator.lt', 'leq' : 'operator.le',
        'ceil' : 'numpy.ceil', 'sqrt' : 'numpy.sqrt',        # libsbml aliases
        'equal' : 'operator.eq', 'not_equal' : 'operator.ne',   # numpy2numpy aliases
        'greater' : 'operator.gt', 'greater_equal' : 'operator.ge', # numpy2numpy aliases
        'less' : 'operator.lt', 'less_equal' : 'operator.le', # numpy2numpy aliases
        'ne' : 'operator.ne', 'ge' : 'operator.ge', 'le' : 'operator.le', # operator2operator
        'xor' : 'operator.xor', 'piecewise' : 'self._piecewise_', '_piecewise_' : 'self._piecewise_',
        'not' : 'operator.not_', 'not_' : 'operator.not_'
    }

    MathmlToNumpy_symb = {
        'notanumber' : 'numpy.NaN', 'pi' : 'numpy.pi',
        'infinity' : 'numpy.Infinity', 'exponentiale' : 'numpy.e',
        'true' : 'True', 'false' : 'False', 'True' : 'True', 'False' : 'False'
    }

    SymbolReplacements = None

    MathmlToInfix = {
        'and' : 'and', 'or' : 'or', 'true' : 'True', 'false' : 'False', 'xor' : 'xor'
    }

    BOOLEANTRUE = ('True','TRUE','true')
    BOOLEANFALSE = ('False','False','false')

    precedence = (
        ('left',  'PLUS', 'MINUS'),
        ('left',  'TIMES', 'DIVIDE'),
        ('left',  'POWER'),
        ('right', 'UMINUS')
        )

    # List of token names
    tokens = ('FIXDEC',
              'FUNCDEC',
              'EVENTDEC',
              'COMPDEC',
              'RRULEDEC',
              'OBJFUNCDEC',
              'FLUXBNDDEC',
              'USERCNSTRDEC',
              'IRREV',
              'REAL',
              'INT',
              'PLUS',
              'MINUS',
              'TIMES',
              'DIVIDE',
              'POWER',
              'LPAREN',
              'RPAREN',
              'EQUALS',
              'SYMBEQUI',
              'COMMA',
              'REACTION_ID',
              'STOICH_COEF',
              'NAME',
              'FORCEFUNC',
              'TIMEFUNC',
              'USERFUNC',
              'INITFUNC',
              'IN',
              'POOL',
              'KEYWORD',
              'KEYWORDB',
              'UNIT',
              'MULTICOMMENT')

    ReservedTerms = RESERVEDTERMS

#################################################################
#################################################################

    def __init__(self,debug=0):
        self.display_debug = debug
        # Lists
        self.LexErrors = []   # List of lexing errors
        self.ParseErrors = []
        self.SymbolErrors = []
        self.ParseOK = True
        self.LexOK = True
        self.ModelUsesNumpyFuncs = False

        self.ReactionIDs = [] # List of reaction names
        self.Names = []       # List of all reagent, parameter and function names
        self.nDict = {}    # Dictionary containing all reaction information
        self.InDict = {}         # spatial dictionary
        self.sDict = {}
        self.pDict = {}    #parameter dict
        self.uDict = {'substance': {'exponent': 1, 'multiplier': 1.0, 'scale': 0, 'kind': 'mole'},
                       'volume': {'exponent': 1, 'multiplier': 1.0, 'scale': 0, 'kind': 'litre'},
                       'time': {'exponent': 1, 'multiplier': 1.0, 'scale': 0, 'kind': 'second'},
                       'length': {'exponent': 1, 'multiplier': 1.0, 'scale': 0, 'kind': 'metre'},
                       'area': {'exponent': 2, 'multiplier': 1.0, 'scale': 0, 'kind': 'metre'}
                      }
        self.KeyWords = {'Modelname' : None,
                         'Description' : None,
                         'Species_In_Conc' : None,
                         'Output_In_Conc' : None
                         }
        self.compartments = {}
        self.InitStrings = []    # Initialisation strings
        self.Inits = []          # Initialised entities
        self.Reagents = []       # All reagents found during parsing of reactions
        self.species = []    # Variable reagents that occur in reactions
        self.FixedReagents = []  # Fixed reagents
        self.ReacParams = []     # Temporary list of reaction parameters
        self.InitParams = []     # Initialised parameters

        self.ForceFunc = []        # Forcing function definition
        self.TimeFunc = []    # Time function definition
        self.UserFunc = []    # User function definition
        self.InitFunc = []    # Intitialization function definition
        self.Functions = {}   # new Function blocks
        self.Events = {}       # new event blocks
        self.AssignmentRules = {}   # new forcing function blocks
        self.UserFuncs = {}   # new user function blocks
        self.ModelInit = {}   # new model initialisation blocks
        self.cbm_FluxBounds = {}  # flux bounds in a CB model
        self.cbm_ObjectiveFunctions = {} # objective functions for a CB model
        self.cbm_UserFluxConstraints = {} # objective functions for a CB model


        self.mach_spec_eps2k = 2.2204460492503131e-14
        self.AllRateEqsGiven = 1 # Flag to check that all rate equations have been given
        self.Debug = 0

        # elementary regular expressions used as building blocks
        self.Int = r'\d+'                                      # Integer
        self.Dec = self.Int + '\.' + self.Int                            # Decimal
        self.Exp = r'([E|e][\+|\-]?)' + self.Int                    # Exponent
        self.Real = self.Dec  + '(' + self.Exp + ')?' + '|' + self.Int + self.Exp  # Real - dec w/o optional exp or int with exp

        # Simple tokens
        self.t_IRREV = r'>'
        self.t_REAL = self.Real
        self.t_INT = self.Int
        self.t_PLUS = r'\+'
        self.t_MINUS = r'-'
        self.t_TIMES = r'\*'
        self.t_DIVIDE = r'/'
        self.t_POWER = '\*\*'
        self.t_LPAREN = r'\('
        self.t_RPAREN = r'\)'
        self.t_EQUALS = r'='
        self.t_COMMA = r','
        self.t_POOL = r'\$pool'
        self.t_IN = r'@'

    t_ignore = ' \t\r'    # Ignore spaces and tabs --- and windows return - brett 20040229

    def t_comment(self,t):
        r'\#.+\n'       # Match from # to newline
        t.lineno += 1   # Increment line number

    def t_multilinecomment(self,t):
        r'"""(.|\n)*?"""'
        t.type = 'MULTICOMMENT'

    def t_newline(self,t):
        r'\n+'          # Match newline
        t.lineno += len(t.value) # Increment with number of consecutive newlines

    def t_SYMBEQUI(self,t):
        r'!=|<'
        t.type = 'SYMBEQUI'
        print('t,{0}'.format(t))
        return t

    def t_FIXDEC(self,t):
        r'FIX:'
        t.type = 'FIXDEC'
        t.value = 'FIX:'
        return t

    def t_KW_Name(self, t):
        r'Modelname:.+\n'
        t.type = 'KEYWORD'
        return t

    def t_KW_Description(self, t):
        r'Description:(.| )+\n'
        t.type = 'KEYWORD'
        return t

    def t_KW_ModelType(self, t):
        r'ModelType:(.| )+\n'
        t.type = 'KEYWORD'
        return t

    def t_KW_Species_In_Conc(self, t):
        r'Species_In_Conc:.+\n'
        t.type = 'KEYWORDB'
        return t

    def t_KW_Output_In_Conc(self, t):
        r'Output_In_Conc:.+\n'
        t.type = 'KEYWORDB'
        return t

    def t_Unit(self, t):
        r'Unit(Substance|Time|Length|Area|Volume):.+\n'
        t.type = 'UNIT'
        return t

    def t_FUNCDEC(self, t):
        r'Function:(.|\n)*?}'
        t.type = 'FUNCDEC'
        return t

    def t_RATERULEDEC(self, t):
        r'RateRule:(.|\n)*?}'
        t.type = 'RRULEDEC'
        return t

    def t_EVENTDEC(self, t):
        r'Event:(.|\n)*?}'
        t.type = 'EVENTDEC'
        return t

    def t_OBJFUNCDEC(self,t):
        r'ObjectiveFunctions:(.|\n)*?}'
        t.type = 'OBJFUNCDEC'
        return t

    def t_FLUXBNDDEC(self,t):
        r'FluxBounds:(.|\n)*?}'
        t.type = 'FLUXBNDDEC'
        return t

    def t_USERCNSTRDEC(self,t):
        r'UserFluxConstraints:(.|\n)*?}'
        t.type = 'USERCNSTRDEC'
        return t


    def t_COMPDEC(self, t):
        r'Compartment:.+\n' # match from cc<space> to end of line
        t.type = 'COMPDEC'
        return t

    def t_FORCEFUNC(self,t):
        r'!F\ .+\n' # match from !F<space> to end of line
        t.type = 'FORCEFUNC'
        t.lineno += 1   # Increment line number
        t.value = t.value[3:]
        return t

    def t_TIMEFUNC(self,t):
        r'!T\ .+\n' # match from !T<space> to end of line
        t.type = 'TIMEFUNC'
        t.lineno += 1   # Increment line number
        t.value = t.value[3:]
        return t

    def t_USERFUNC(self,t):
        r'!U\ .+\n' # match from !U<space> to end of line
        t.type = 'USERFUNC'
        t.lineno += 1   # Increment line number
        t.value = t.value[3:]
        return t

    def t_INITFUNC(self,t):
        r'!I\ .+\n' # match from !I<space> to end of line
        t.type = 'INITFUNC'
        t.lineno += 1   # Increment line number
        t.value = t.value[3:]
        return t

    def t_REACTION_ID(self,t):
        r'[a-zA-Z_][\w@]*:' # Match any letter in [a-zA-Z0-9_] up to a colon
                            # allow "_" to startswith - brett 20050823
        t.type = 'REACTION_ID'
        t.value = t.value[:-1]  # remove the colon

        #print t.value, self.ReactionIDs
        #sleep(1)
        if '@' in t.value:
            ts = t.value.split('@')
            t.value = ts[0]
            self.InDict.update({ts[0] : ts[1]})
        if t.value in self.ReactionIDs:
            # brett I think this is hyperactive reporting ... removing 20050321
            # self.LexErrors.append(('Duplicate ReactionID ', t.lineno, t.value, t.type))
            pass
        else:
            self.ReactionIDs.append(t.value)
        # possible alternative to above - brett
        #if t.value not in self.ReactionIDs:
        #    self.ReactionIDs.append(t.value)
        return t

    def t_STOICH_COEF(self,t):
        r'\{\d+\}|\{\d+\.\d+\}'
        t.type = 'STOICH_COEF'
        t.value = t.value[1:-1]  # Remove left and right curly brackets
        return t

    def t_NAME(self,t):
        r'numpy\.[\w]*|math\.[\w]*|operator\.[\w]*|[a-zA-Z_][\w@]*'
        SciCons = False
        if '@' in t.value:
            ts = t.value.split('@')
            t.value = ts[0]
            self.InDict.update({ts[0] : ts[1]})
        if t.value in self.MathmlToNumpy_symb:
            if self.MathmlToNumpy_symb[t.value] == None:
                self.SymbolErrors.append(t.value)
                print('\nSymbol \"{0}\" not yet supported by PySCeS.'.format(t.value))
                gt = 'unknown_symbol_' + t.value
                t.value = gt
            else:
                gt = self.MathmlToNumpy_symb[t.value]
                t.value = gt
            self.ModelUsesNumpyFuncs = 1
            SciCons = True
        elif self.SymbolReplacements != None and t.value in self.SymbolReplacements:
            if t.value not in self.Names:
                self.Names.append('self.' + self.SymbolReplacements[t.value])
            gt = 'self.' +  self.SymbolReplacements[t.value]
            t.value = gt
        elif t.value not in self.Names and t.value != '_TIME_': # Only add to list if absent in list
            self.Names.append('self.' + t.value)
        if t.value[:6] == 'numpy.' or t.value[:5] == 'math.' or t.value[:9] == 'operator.':
            pass
        elif t.value not in self.MathmlToNumpy_funcs and not SciCons: # make class attributes, ignore function names
            gt = 'self.' + t.value
            t.value = gt
        t.type = 'NAME'
        return t

    def t_error(self,t):
        ##  try:
            ##  self.LexErrors.append(('Lexer error ', t.lineno, t.value, t.type))
        ##  except:
            ##  print 'Lexer error'
        ##  #print 'Illegal character, Line ' + str(t.lineno) + ' :' + str(t.value[0])
        ##  t.skip(1)

        print("Illegal character '{0}'".format(t.value[0]))
        self.LexErrors.append(t.value[0])
        self.LexOK = False
        t.lexer.skip(1)

    ##############
    # The parser #
    ##############

    def Show(self,name,tok):
        if self.Debug:
            print(name,tok)

    def p_error(self,t):
        try:
            ##  self.ParseErrors.append(('Syntax error ', t.lineno, t.value, t.type))
            self.ParseErrors.append(t)
        except:
            print('p_error generated a parsing error')
        tok = yacc.token()
        while tok and tok.type != 'REACTION_ID':
            tok = yacc.token()
        self.ParseOK = False
        return tok

    def p_model(self,t):
        '''Model : Statement
                 | Model Statement '''

        self.Show('Model',t[0])

    def p_statement(self,t):
        '''Statement : Fixed
                     | FunctionDec
                     | RateRuleDec
                     | EventDec
                     | ObjFuncDec
                     | FluxBndDec
                     | UserCnstrDec
                     | CompartmentDec
                     | ReactionLine
                     | Initialise
                     | Forcedfunc
                     | Timefunc
                     | Userfunc
                     | Initfunc
                     | NameInName
                     | KeyWord
                     | KeyWordB
                     | Unit
                     | MultiComment
                     | SymbEqui'''
        self.Show('Statement',t[0])

    def p_nameinname(self, t):
        '''NameInName : NAME IN NAME'''
        print('This should never fire, what we need to do is ...')

    def p_inequalities_symb(self, t):
        '''SymbEqui : SYMBEQUI'''
        print('p_inequalities_symb', t[:])
        t[0] = t[1]

    def p_multicomment(self, t):
        '''MultiComment : MULTICOMMENT'''
        self.Show('KeyWord:',t[0])

    def p_keyword(self, t):
        '''KeyWord : KEYWORD'''
        ##  print 'KeyWord:', t[:]
        kpair = [e.strip() for e in t[1].split(':')]
        if kpair[0] in self.KeyWords:            
            self.KeyWords.update({kpair[0] : kpair[1]})
        ##  print self.KeyWords
        self.Show('KeyWord:',t[0])

    def p_keywordboolean(self, t):
        '''KeyWordB : KEYWORDB'''
        ##  print 'KeyWordB:', t[:]
        kpair = [e.strip() for e in t[1].split(':')]        
        if kpair[0] in self.KeyWords:            
            if kpair[1] in self.BOOLEANTRUE:
                self.KeyWords.update({kpair[0] : True})
            elif kpair[1] in self.BOOLEANFALSE:
                self.KeyWords.update({kpair[0] : False})
        ##  print self.KeyWords
        self.Show('KeyWordB:',t[0])

    def p_unit(self, t):
        '''Unit : UNIT'''
        ##  print 'u', t[:]
        u = t[1].split(',')
        u[0] = u[0].split(':')
        u.append(u[0][1])
        u[0] = u[0][0].lower()
        ##  print u
        u[0] = u[0].replace('unit','')
        ##  print u
        for i in range(len(u)):
            u[i] = u[i].strip()
            if i in [1,2,3]:
                u[i] = float(u[i])
        ##  print u
        self.uDict.update({u[0] : {'multiplier': u[1],
                                    'scale': u[2],
                                    'exponent': u[3],
                                    'kind': u[4]
                                    }
                           })
        self.Show('Unit:',t[0])

    def p_fixed(self,t):
        '''Fixed : FIXDEC FixList'''
        self.Show('Fixed:',t[0])

    def p_functiondec(self, t):
        '''FunctionDec : FUNCDEC'''
        ##  print 'p_functiondec', t[:]
        rawf = t[1].replace('Function:','').lstrip()
        args = rawf[:rawf.find('{')].strip().split(',')
        name = args.pop(0)
        func = rawf[rawf.find('{')+1:rawf.find('}')]
        self.Functions.update({name : {
                                      'name' : name,
                                      'args' : args,
                                      'formula' : func
                                    }
                            })
        self.Show('FunctionDec:',t[0])

    def p_rateruledec(self, t):
        '''RateRuleDec : RRULEDEC'''
        ##  print 'p_functiondec', t[:]
        rawf = t[1].replace('RateRule:','').lstrip()
        name = rawf[:rawf.find('{')].strip()
        func = rawf[rawf.find('{')+1:rawf.find('}')]
        self.AssignmentRules.update({name : {'name' : name,
                                              'formula' : func.strip(),
                                              'type' : 'rate'
                                            }
                                    })
        self.Show('RateRuleDec:',t[0])

    def p_eventdec(self, t):
        '''EventDec : EVENTDEC'''
        rawf = t[1].replace('Event:','').lstrip()
        args = rawf[:rawf.find('{')].strip().split(',')
        name = args.pop(0)
        delay = float(args.pop(-1))
        trigger = ''
        for a in args:
            trigger = trigger + a + ','
        trigger = trigger[:-1]

        rawF = rawf[rawf.find('{')+1:rawf.find('}')].split('\n')
        assignments = {}
        for ass in rawF:
            if len(ass.strip()) > 0:
                ass = ass.split('=')
                assignments.update({ass[0].strip() : ass[1].strip()})
        self.Events.update({name : { 'delay' : delay,
                                      'name' : name,
                                      'trigger' : trigger,
                                      'assignments' : assignments,
                                      'tsymb' : None
                                    }
                            })
        self.Show('EventDec:',t[0])

    def p_objfuncdec(self, t):
        '''ObjFuncDec : OBJFUNCDEC'''
        rawf = t[1].replace('ObjectiveFunctions:','').lstrip()
        args = rawf[:rawf.find('{')].strip().split(',')
        rawf = [r.strip() for r in rawf[rawf.find('{')+1:rawf.find('}')].split('\n')]
        objectives = {}
        for oo in rawf:
            print(rawf)
            print(oo)
            if len(oo.strip()) > 0:
                oS = [t.strip() for t in oo.split(':')]
                active = False
                if oS[0] == args[0]:
                    active = True
                objectives.update({oS[0] : {'id' : oS[0],
                                             'active' : active,
                                             'osense' : oS[1],
                                             'fluxes' : oS[2]}
                                            })
        self.cbm_ObjectiveFunctions.update(objectives)
        self.Show('ObjFuncDec:',t[0])

    def p_fluxbnddec(self, t):
        '''FluxBndDec : FLUXBNDDEC'''
        rawf = t[1].replace('FluxBounds:','').lstrip()
        rawf = [r.strip() for r in rawf[rawf.find('{')+1:rawf.find('}')].split('\n')]

        fluxbnds = {}
        cntr = 0
        for oo in rawf:
            if len(oo.strip()) > 0:
                operator = None
                lhs = None
                rhs = None
                fid = 'fbnd_{0}'.format(cntr)
                for oper in ['>=','<=','>','<','=']:
                    if  oper in oo:
                        bnd = oo.split(oper)
                        operator = oper
                        lhs = bnd[0].strip()
                        rhs = float(bnd[1].strip())
                        break
                fluxbnds.update({fid : {'id' : fid,
                                           'flux' : lhs,
                                           'rhs' : rhs,
                                           'operator' : operator}
                                          })
                cntr += 1
        self.cbm_FluxBounds.update(fluxbnds)
        self.Show('FluxBndDec:',t[0])

    def p_usercnstrdec(self, t):
        '''UserCnstrDec : USERCNSTRDEC'''
        rawf = t[1].replace('UserFluxConstraints:','').lstrip()
        rawf = [r.strip() for r in rawf[rawf.find('{')+1:rawf.find('}')].split('\n')]

        fluxcnstr = {}
        for oo in rawf:
            if len(oo.strip()) > 0:
                oS = [t.strip() for t in oo.split(':')]
                id = oS[0]
                operator = None
                lhs = None
                rhs = None
                for oper in ['>=','<=','>','<','=']:
                    if  oper in oS[1]:
                        cnstr = oS[1].split(oper)
                        operator = oper
                        lhs = [c.strip() for c in cnstr[0].split(',')]
                        lhs2 = []
                        for c2 in lhs:
                            c = c2.split(' ')
                            if len(c) == 1:
                                coeff = '+1'
                                val = c[0]
                            else:
                                coeff = c[0]
                                val = c[1]
                            lhs2.append((float(coeff),val))
                        rhs = float(cnstr[1].strip())
                        break
                fluxcnstr.update({oS[0] : {'id' : oS[0],
                                             'constraint' : tuple(lhs2),
                                             'rhs' : rhs,
                                             'operator' : operator}
                                            })
        self.cbm_UserFluxConstraints.update(fluxcnstr)
        self.Show('UserCnstrDec:',t[0])




        self.Show('UserCnstrDec:',t[0])





    def p_compartmentdec(self,t):
        '''CompartmentDec : COMPDEC'''
        rawf = t[1].strip().replace('Compartment:','')
        strAr = [st.strip() for st in rawf.split(',')]
        if len(strAr) <= 3:
            name = strAr[0]
            size = strAr[1]
            dimensions = strAr[2]
        if len(strAr) == 4:
            area = strAr[3]
        else:
            area = None

        self.compartments.update({name:{'name':name,
                                         'size': float(size),
                                         'dimensions' : int(dimensions),
                                         'compartment': None,
                                         'area' : None
                                         ##  'volume' : None
                                          }})

    def p_forcedfunc(self,t):
        '''Forcedfunc : FORCEFUNC'''
        self.ForceFunc.append(t[1])
        ar = t[1].split('=')
        name = ar[0].replace('self.','').strip()
        func = ar[1].replace('self.','').strip()
        self.AssignmentRules.update({name : {'name' : name,
                                              'formula' : func,
                                              'type' : 'assignment'
                                            }
                                    })
        self.Show('Forcedfunc:',t[0])

    def p_timefunc(self,t):
        '''Timefunc : TIMEFUNC'''
        self.TimeFunc.append(t[1])
        self.Show('Timefunc:',t[0])

    def p_userfunc(self,t):
        '''Userfunc : USERFUNC'''
        self.UserFunc.append(t[1])
        ar = t[1].split('=')
        name = ar[0].replace('self.','').strip()
        func = ar[1].replace('self.','').strip()
        self.UserFuncs.update({name : {'name' : name,
                                            'formula' : func,
                                            'type' : 'userfuncs'}
                                    })
        self.Show('Userfunc:',t[0])

    def p_initfunc(self,t):
        '''Initfunc : INITFUNC'''
        self.InitFunc.append(t[1])
        ar = t[1].split('=')
        name = ar[0].replace('self.','').strip()
        value = ar[1].replace('self.','').replace('\'','').replace('"','').strip()
        self.ModelInit.update({name : value})
        self.Show('Initfunc:',t[0])

    def p_fixedreagents(self,t):
        '''FixList : NAME
                   | NAME FixList'''
        if t[1] != None:
            self.FixedReagents.append(t[1])
        t[0] = [t[1]]
        try:
            t[0] += t[2]
        except:
            pass
        self.Show('FixList', t[0])

    # TODO:
    def p_initialise(self,t):
        '''Initialise : NAME EQUALS Expression'''
        value = None
        name = t[1].replace('self.','')
        try:
            value = eval(t[3])
            # 20090402 we need to keep track of species parameter initialisations and then
            # perform a lookup or implement a proxy class that can store
            # evaluated assignments
        except Exception as ex:
            print('Initialisation error: {0}'.format(t[3]) )
            print(ex)

        self.sDict.update({name : {'name' : name,
                                    'initial' : value
                                   }
                           })
        t[0] = t[1] + t[2] + t[3]

        self.InitStrings.append(t[0])
        self.Inits.append(t[1])
        self.Show('Initialisation',t[0])

    def p_reaction_line(self,t):
        '''ReactionLine : REACTION_ID ReactionEq
                        | REACTION_ID ReactionEq Expression'''

        #global AllRateEqsGiven, self.ReacParams
        ReacID = t[1]
        if ReacID in self.nDict:
            self.ParseErrors.append(('Duplicate Reaction ', t.lineno, ReacID, None))
        self.nDict[ReacID] = {} # Reaction dictionary for ReacID
        ##  self.nDict[ReacID]['Reagents'] = {} # Reagent dictionary within ReacID
        self.nDict[ReacID].update({'Reagents' : {}, 'name' : ReacID})
        # a list of all reagents specified in the stoichiometry, unmodified and not taken into consideration for anything
        self.nDict[ReacID].update({'AllReagents' : t[2][0]})

        # brett: if an index exists sum the coefficients instead of adding a new one
        # this seems to deal with multiple definitions like X + X > Y and 2{X} + Y > Z + X
        for i in t[2][0]: # First tuple member of ReactionEq contains list of (name,stoichcoef)
            if i[0] in self.nDict[ReacID]['Reagents']:
                self.nDict[ReacID]['Reagents'][i[0]] = self.nDict[ReacID]['Reagents'][i[0]] + i[1]
            else:
                self.nDict[ReacID]['Reagents'][i[0]] = i[1] # Key for reagent with stoichcoef value
        killList = []

        # brett: however for the case of X + Y > Y + Z where the sum of the coefficients
        # is zero we can delete the key (Y) out of the reaction list altgether (hopefully!)
        for i in self.nDict[ReacID]['Reagents']:
            if abs(self.nDict[ReacID]['Reagents'][i]) < self.mach_spec_eps2k:
                killList.append(i)
                #print self.mach_spec_eps2k, self.nDict[ReacID]['Reagents']
        #print killList, self.nDict[ReacID]['Reagents']
        # brett: and the easiest way of doing this is putting the zero keys in a list
        # and deleting them out of the dictionary
        if len(killList) != 0:
            for i in killList:
                del self.nDict[ReacID]['Reagents'][i]
        #print killList, self.nDict[ReacID]['Reagents']

        self.nDict[ReacID]['Type'] = t[2][1] # Second tuple member of ReactionEq contains type
        try: # Save rate equation and create parameter list
            self.nDict[ReacID]['RateEq']   = t[3]
            self.nDict[ReacID]['Params']   = self.ReacParams
            self.ReacParams = [] # Reset global self.ReacParams list
        except:
            self.nDict[ReacID]['RateEq']   = ''
            self.nDict[ReacID]['Params']   = []
            self.AllRateEqsGiven = 0 # Set global flag to false
        self.Show('ReactionLine',t[0])

    def p_reaction_eq(self,t):
        '''ReactionEq : LeftHalfReaction EQUALS RightHalfReaction
                      | LeftHalfReaction IRREV  RightHalfReaction
                      | POOL EQUALS  RightHalfReaction
                      | POOL IRREV  RightHalfReaction
                      | LeftHalfReaction EQUALS POOL
                      | LeftHalfReaction IRREV POOL'''

        ReacType = ''
        if   t[2] == '=':
            ReacType = 'Rever'
        elif t[2] == '>':
            ReacType = 'Irrev'

        # Yeah well you know the story ... 4 hrs of trying the "cool, other" way of doing it at work
        # and then solving the original problem in 90 mins with 14 lines of code at home ...
        # oh almost forgot: anonymous source and sink pools now work in PySCeS - brett 20050908
        if t[1] == '$pool':
            t[0] = (t[3], ReacType)
        elif t[3] == '$pool':
            t[0] = (t[1], ReacType)
        else:
            t[0] = (t[1] + t[3], ReacType)
        #print 'reaction_eq',t[0]
        self.Show('ReactionEq',t[0])

    def p_left_half_reaction(self,t):
        ''' LeftHalfReaction : SubstrateTerm
                             | SubstrateTerm PLUS LeftHalfReaction'''
        # Make a list of substrate terms
        t[0] = [t[1]]
        try:
            t[0] += t[3]
        except:
            pass
        self.Show('LeftHalfReaction', t[0])

    def p_right_half_reaction(self,t):
        ''' RightHalfReaction : ProductTerm
                              | ProductTerm PLUS RightHalfReaction'''

        t[0] = [t[1]]

        try:
            t[0] += t[3]
        except:
            pass
        self.Show('RightHalfReaction', t[0])

    def p_substrate_term(self,t):
        '''SubstrateTerm : STOICH_COEF NAME
                         | NAME'''

        # Make tuple of NAME and stoichiometric coefficient
        # (< 0 because substrate)
        try:
            t[0] = (t[2], -float(t[1]))
            if t[2] not in self.Reagents:# and t[2] != 'self.pool':
                self.Reagents.append(t[2])
        except:
            t[0] = (t[1], -1.0)
            if t[1] not in self.Reagents:#  and t[1] != 'self.pool':
                self.Reagents.append(t[1])
        self.Show('SubstrateTerm', t[0])

    def p_product_term(self,t):
        '''ProductTerm : STOICH_COEF NAME
                       | NAME'''
        # Make tuple of NAME and stoichiometric coefficient
        # (> 0 because product)
        try:
            t[0] = (t[2], float(t[1]))
            if t[2] not in self.Reagents:# and t[2] != 'self.pool':
                self.Reagents.append(t[2])
        except:
            t[0] = (t[1], 1.0)
            if t[1] not in self.Reagents:# and t[1] != 'self.pool':
                self.Reagents.append(t[1])
        self.Show('ProductTerm', t[0])

    def p_rate_eq(self,t):
        '''Expression : Expression PLUS Expression
                      | Expression MINUS Expression
                      | Expression TIMES Expression
                      | Expression DIVIDE Expression
                      | Power
                      | Number
                      | Func'''
                    # |UMINUS : add if the
                    #  alternative for p_uminus is used

        if len(t.slice)==4:
            t[0] = t[1] + t[2] + t[3]
        else:
            t[0] = t[1]

    def p_power(self,t):
        '''Power : Expression POWER Expression'''

        t[0] = 'pow('+ t[1] + ',' + t[3] + ')' #changed to make it DeriVar compatible
        ##  t[0] = t[1] + '**' + t[3]

    def p_uminus(self,t):
        '''Expression : MINUS Expression %prec UMINUS'''
        # Alternative '''UMINUS : MINUS Expression'''

        t[0] = t[1] + t[2]

    def p_number(self,t):
        '''Number : REAL
                  | INT
                  | NAME'''

    # Build list of entities
        ttype = 'param'
        try:
            tx = float(t[1])
            ttype = 'float'
        except ValueError:
            pass
            ##  try:
                ##  t[1] = complex(t[1])
                ##  ttype = 'complex'
            ##  except:
                ##  pass
        if ttype in ['complex','float']:
            t[1] = str(tx)
        else:
            # if this is not a number add it as a parameter EXCEPT if it is a function's
            # name or numpy.constant
            if (t[1] not in self.ReacParams) and\
                (t[1].replace('numpy.','').replace('math.','').replace('operator.','') not in self.MathmlToNumpy_funcs) and\
                (t[1].replace('numpy.','').replace('math.','').replace('operator.','') not in self.MathmlToNumpy_symb): # ignore function names and duplications
                self.ReacParams.append(t[1])
        t[0] = t[1]

     # new Core2 based InfixParser code
    def p_function(self,t):
        '''Func : LPAREN ArgList RPAREN
                | NAME LPAREN ArgList RPAREN'''

        # convert root(degree,<expr>) to pow(<expr>, 1/degree)
        if t[1].strip() == 'root':
            t[1] = self.MathmlToNumpy_funcs[t[1]]
            t[3] = '{0}, {1}'.format(t[3][t[3].index(',')+1:], 1.0/float(t[3][:t[3].index(',')])  )
            t[0] = t[1] + t[2] + t[3] + t[4]
        elif t[1] in self.MathmlToNumpy_funcs:
            if self.MathmlToNumpy_funcs[t[1]] == None:
                self.SymbolErrors.append(t[1])
                print('\nFunction \"{0}\" not supported'.format(t[1]))
                t[0] = 'unknown_function_'+t[1] + t[2] + t[3] + t[4]
            else:
                try:
                    t[0] = self.MathmlToNumpy_funcs[t[1]] + t[2] + t[3] + t[4]
                except Exception as EX:
                    print('Function Parse error 1 (please report!)\n', EX)
            self.ModelUsesNumpyFuncs = True
        # differentiate between bracketed functions and expressions:
        # func( S1 ) and ( S/S05 )
        elif len(t) == 4:
            t[0] = t[1] + t[2] + t[3]
        else:
            # assume some arbitrary function definition
            if t[1][:6] == 'numpy.' or t[1][:5] == 'math.' or t[1][:9] == 'operator.': # NEW UNTESTED
                t[0] = t[1] + t[2] + t[3] + t[4]
            else:
                t[0] = t[1] + t[2] + t[3] + t[4]

    # arbitrary number of arguments in an argument list
    # adapted from Andrew Dalke's GardenSnake
    # http://www.dalkescientific.com/writings/diary/GardenSnake.py
    def p_arglist(self,t):
        '''ArgList : Expression
                   | ArgList COMMA Expression'''
        try:
            if len(t) == 2:
                t[0] = t[1]
            elif len(t) == 4:
                t[0] = t[1] + ',' + t[3]
            else:
                pass
        except Exception as EX:
            print('Function ArgList error (please report!)\n', EX)

###################################################################
###################################################################

    def ParsePSC(self,modelfile,modeldir,modeloutput,quiet):
        """
        ParsePSC(modelfile,modeldir,modeloutput)

        Parse the .psc file into a Network Dictionary and associated property arrays

        modelfile: PSC filename
        modeldir: PSC input directory
        modeloutput: working directory for lex/parse temporary files
        quiet: quite mode

        """
        # we clear stuff so we have a shiny new instance
        self.nDict = {}    # Dictionary containing all reaction information
        self.InDict = {}
        self.sDict = {}
        self.pDict = {}    #parameter dict
        self.uDict = {'substance': {'exponent': 1, 'multiplier': 1.0, 'scale': 0, 'kind': 'mole'},
                       'volume': {'exponent': 1, 'multiplier': 1.0, 'scale': 0, 'kind': 'litre'},
                       'time': {'exponent': 1, 'multiplier': 1.0, 'scale': 0, 'kind': 'second'},
                       'length': {'exponent': 1, 'multiplier': 1.0, 'scale': 0, 'kind': 'metre'},
                       'area': {'exponent': 2, 'multiplier': 1.0, 'scale': 0, 'kind': 'metre'}
                      }
        self.KeyWords = {'Modelname' : None,
                         'Description' : None,
                         'Species_In_Conc' : None,
                         'Output_In_Conc' : None,
                         'ModelType' : None
                         }
        self.compartments = {}
        self.ReactionIDs = [] # List of reaction names
        self.InitStrings = []    # Initialisation strings
        self.species = []    # Variable reagents that occur in reactions
        self.FixedReagents = []  # Fixed reagents
        self.InitParams = []     # Initialised parameters
        self.Names = []       # List of all reagent, parameter and function names
        self.Inits = []          # Initialised entities
        self.Reagents = []       # All reagents found during parsing of reactions
        self.ReacParams = []     # Temporary list of reaction parameters
        self.LexErrors = []   # List of lexing errors
        self.ParseErrors = []
        self.SymbolErrors = []
        self.ForceFunc = []
        self.TimeFunc = []
        self.UserFunc = []
        self.InitFunc = []
        self.Functions = {}
        self.Events = {}
        self.AssignmentRules = {}
        self.UserFuncs = {}
        self.ModelInit = {}
        self.cbm_FluxBounds = {}  # flux bounds in a CB model
        self.cbm_ObjectiveFunctions = {} # objective functions for a CB model
        self.cbm_UserFluxConstraints = {} # objective functions for a CB model

        self.ModelUsesNumpyFuncs = False

        # 4 hrs of debugging and these two lines solve the problem .... I'm out of here!
        importlib.reload(lex)
        importlib.reload(yacc)
        # fixes the obscure reload <--> garbage collection reference overload bug ... who said scripted lang's were
        # easy to use :-) - brett 20040122

        self.AllRateEqsGiven = 1
        self.ParseOK = True
        self.LexOK = True

        # check to see if the last line of the file is an os defined empty line
        self.CheckLastLine(os.path.join(modeldir,modelfile))

        # load model data from file
        Data = open(os.path.join(modeldir,modelfile),'r')
        self.__Model = Data.read()
        Data.close()

        # try and find a temporary workspace or use ModelOutput
        #if 'TMP' in os.environ:
        #    tempDir = os.environ['TMP']
        #elif 'TEMP' in os.environ:
        #    tempDir = os.environ['TEMP']
        #elif os.path.isdir(modeloutput):
        #    tempDir = modeloutput
        #else:
        #    tempDir = os.getcwd()

        os.chdir(output_dir) # tempDir before

        # fix filenames for intermediary files - brett
        if not modelfile[:-4].isalnum():
            FileL = list(modelfile)
            FileT = ''
            for let in FileL:
                if let.isalnum():
                    FileT += let
            self.debugfile = '_pys' + FileT[:-3] + ".dbg"
            #self.tabmodule = '_pys' + FileT[:-3] + "_" + "parsetab"
            self.tabmodule = '_pys'  "_" + "parsetab"
        else:
            self.debugfile = '_pys' + modelfile[:-4] + ".dbg"
            #self.tabmodule = '_pys' + modelfile[:-4] + "_" + "parsetab"
            self.tabmodule = '_pys' + "_" + "parsetab"

        if self.Debug:
            print(self.tabmodule)
            print(self.debugfile)

        self.__lexer = lex.lex(module=self, debug=self.Debug)
        self.__lexer.input(self.__Model)
        self.__parser = yacc.yacc(module=self,
                debug=self.Debug,
                debugfile=self.debugfile,
                tabmodule=self.tabmodule)

        while 1:
            tok = self.__lexer.token()
            if not tok: break

        while 1:
            p = self.__parser.parse(self.__Model)
            if not p: break

        # get to our work directory
        os.chdir(modeloutput)

        # we have the dictionary get rid of this stuff
        del self.__Model, self.__lexer, self.__parser, p

        # Create forcing function code blocks - brett 20050621
        Fstr = ''
        Tstr = ''
        Ustr = ''
        Istr = ''
        for f in self.ForceFunc:
            if os.sys.platform != 'win32' and f[-2] == '\r':
                Fstr = Fstr + f[:-2] + '\n'
            else:
                Fstr += f
        for t in self.TimeFunc:
            if os.sys.platform != 'win32' and t[-2] == '\r':
                Tstr = Tstr + t[:-2] + '\n'
            else:
                Tstr += t
        for u in self.UserFunc:
            if os.sys.platform != 'win32' and u[-2] == '\r':
                Ustr = Ustr + u[:-2] + '\n'
            else:
                Ustr += u
        for i in self.InitFunc:
            if os.sys.platform != 'win32' and i[-2] == '\r':
                Istr = Istr + i[:-2] + '\n'
            else:
                Istr += i

        self.ForceFunc = Fstr
        self.TimeFunc = Tstr
        self.UserFunc = Ustr
        self.InitFunc = Istr
        del Fstr,Tstr,Ustr,Istr

        # Create list of variable species
        for reag in self.Reagents:
            if reag not in self.FixedReagents:
                self.species.append(reag)

        # Warn if no reagents have been fixed
        if not quiet:
            if self.FixedReagents == []:
                print('Info: No reagents have been fixed')
            else: # Warn if a fixed reagent does not occur in a reaction equation
                for reag in self.FixedReagents:
                    if reag not in self.Reagents:
                        print('Warning: species "{0:s}" (fixed) does not occur in any reaction'.format(reag.replace('self.','')) )

        # for the initialised ones
        for s in list(self.sDict):
            if 'self.'+s in self.FixedReagents+self.Reagents:
                fixed = False
                compartment = None
                isamount = False
                initial = self.sDict[s]['initial']
                name = self.sDict[s]['name']

                if 'self.'+s in self.FixedReagents:
                    fixed = True
                if s in list(self.InDict):
                    compartment = self.InDict[s]
                self.sDict.update({s : {'name' : name,
                                        'initial' : initial,
                                        'compartment' : compartment,
                                        'fixed' : fixed,
                                        'isamount' : isamount
                                        }})
            else:
                self.pDict.update({s : self.sDict.pop(s)})


        # for the uninitialised ones
        ##  print self.FixedReagents, self.Reagents, self.Inits
        for s in self.species+self.FixedReagents:
            if not s.replace('self.','') in self.sDict:
                fixed = False
                compartment = None
                isamount = False
                initial = None
                name = s.replace('self.','')

                if s in self.FixedReagents:
                    fixed = True
                if name in list(self.InDict):
                    compartment = self.InDict[name]
                self.sDict.update({name : {'name' : name,
                                        'initial' : initial,
                                        'compartment' : compartment,
                                        'fixed' : fixed,
                                        'isamount' : isamount
                                        }})


        # Create list of initialised parameters
        # eventually this must be switched over to pDict
        for elem in self.Inits:
            el = elem.replace('self.', '')
            names = [self.sDict[s]['name'] for s in self.sDict]
            if el not in names:
                self.InitParams.append(elem)
            if el in names and self.sDict[el]['fixed'] == True:
                self.InitParams.append(elem)

        # In self.nDict, clean rate equation pameter list of variables that occur in that reaction
        for id in list(self.nDict):
            # create 'Modifiers' key for each reaction - brett 20050606
            self.nDict[id].update({'Modifiers':[]})
            reagcomp = None
            if id in list(self.InDict):
                reagcomp = self.InDict[id]
            self.nDict[id].update({'compartment' : reagcomp})
            for reag in self.species:
                if reag in self.nDict[id]['Params']:
                    # if the reagent is a reaction reagent delete - brett 20050606
                    if reag in self.nDict[id]['Reagents']:
                        self.nDict[id]['Params'].remove(reag)
                    # if not it is a modifier ... add to modifiers and delete - brett 20050606
                    else:
                        #print 'INFO: variable modifier \"' + reag.replace('self.','') + '\" detected in', id
                        self.nDict[id]['Modifiers'].append(reag)
                        self.nDict[id]['Params'].remove(reag)

        # Check whether all parameters have been initialised and if not add to dictionary
        for id in list(self.nDict):
            cnames = [self.compartments[s]['name'] for s in self.compartments]
            for param in self.nDict[id]['Params']:
                ##  print param
                ##  print self.InitParams
                if param not in self.InitParams and param.replace('self.','') not in cnames:
                    ##  print param, self.InitParams, self.pDict
                    print('Warning: {0:s} parameter "{1:s}"  has not been initialised'.format (id, param.replace('self.','')))
                    pname = param.replace('self.','')
                    self.pDict.update({pname : {'name':pname, 'initial':None}})
                    self.InitParams.append(param)

        # Check whether all variable reagents have been initialised
        for reag in self.species+self.FixedReagents:
            if reag not in self.Inits and reag in self.species:
                print('Warning: species "{0:s}" has not been initialised'.format(reag.replace('self.','')))
            if reag not in self.Inits and reag in self.FixedReagents:
                print('Warning: species "{0:s}" (fixed) has not been initialised'.format(reag.replace('self.','')))

        # Check that all initialised parameters actually occur in self.Inits
        known = 0
        for param in self.InitParams:
            for id in list(self.nDict):
                if param in self.nDict[id]['Params']:
                    known = 1
                    break
                else:
                    known = 0
            if not known: print('Info: "{0:s}" has been initialised but does not occur in a rate equation'.format(param.replace('self.','')))

        # find modifiers for each reaction - brett 20050606
        self.modifiers = []
        for i in self.nDict:
            for j in self.nDict[i]['Params']:
                if j in self.FixedReagents and j not in self.nDict[i]['Reagents']:
                    self.nDict[i]['Modifiers'].append(j)
            self.modifiers.append((i,tuple([j.replace('self.','') for j in self.nDict[i]['Modifiers']])))

        self.fixed_species = copy.copy(self.FixedReagents)
        self.parameters = copy.copy(self.InitParams)
        self.reactions = copy.copy(self.ReactionIDs)

        for x in range(len(self.fixed_species)):
            self.fixed_species[x] = self.fixed_species[x].replace('self.','')
        for x in range(len(self.species)):
            self.species[x] = self.species[x].replace('self.','')
        for x in range(len(self.parameters)):
            self.parameters[x] = self.parameters[x].replace('self.','')

        if self.display_debug == 1:
            print(self.fixed_species)
            print(self.species)
            print(self.parameters)
            print(self.reactions)

        if self.display_debug == 1:
            print('\nDebugging Part 1')
            print('\nFixedReagents')
            print(self.FixedReagents)
            print('\nInitStrings')
            print(self.InitStrings)
            print('\nInitParams')
            print(self.InitParams)
            print('\nReactionIDs')
            print(self.ReactionIDs)
            print('\nspecies')
            print(self.species)
            print('\nNetworkDict')
            print(self.nDict)
            print('\n\n')

    def KeywordCheck(self,inputarr,bad=[]):
        """
        KeywordCheck(inputarr,bad=[])

        Check a list of names for reserved PySCeS keywords

        Arguments:
        =========
        inputarr: a list of names to check
        bad [default=[]]: a list of illegal keywords (new if not supplied

        """
        for x in inputarr:
            if x.replace('self.','') in self.ReservedTerms:
                bad.append(x.replace('self.',''))
                print(x.replace('self.',''), '\t\tERROR: keyword')
        return bad

    def CheckLastLine(self,name):
        """
        CheckLastLine(name)

        Checks to see if file ends with an empty line ... if not it adds one

        Arguments:
        =========
        name: filename of the PySCeS input file

        """
        
        try: F = open(name,'rb+')
        except:
           raise IOError("*** File does not exist ***") # updated 05 08 2014 (read-in as binary mode) for python 3 comp.
        
        F.seek(-1,2)        
        if F.read().decode(encoding='UTF-8') not in ['\n']:
            if os.sys.platform == 'win32':
                F.read()
            try:
                F.write('\n'.encode('utf-8'))
            except Exception as EX:
                print('EOL add error', EX)
            print('Adding \\n to input file')
        F.close()
##      # you REALLY do NOT want to know what this was for! - brett
##      if gc.isenabled() == 1:
##          gc.collect()
##          del gc.garbage[:]
##          gc.set_debug(0)
##          print 'gc debug level = ' + `gc.get_debug()`

