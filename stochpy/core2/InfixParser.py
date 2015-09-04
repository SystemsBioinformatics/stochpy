"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2012 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa.

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license. See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""
from __future__ import print_function
from __future__ import absolute_import

import os,math,imp

from .version import __version__
from . import  lex
from . import yacc

class MyInfixLexer:
    """
    This parser has been cranked to handle Python infix, numpy and MathML2.0 prefix expressions
    """
    debug = 0
    LexOK = True
    LexErrors = None
    __pwcntr__ = 0

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
        'ceil' : 'numpy.ceil', 'sqrt' : 'math.sqrt',        # libsbml aliases
        'equal' : 'operator.eq', 'not_equal' : 'operator.ne',   # numpy2numpy aliases
        'greater' : 'operator.gt', 'greater_equal' : 'operator.ge', # numpy2numpy aliases
        'less' : 'operator.lt', 'less_equal' : 'operator.le', # numpy2numpy aliases
        'ne' : 'operator.ne', 'ge' : 'operator.ge', 'le' : 'operator.le', # operator2operator
        'piecewise' : 'self._piecewise_', '_piecewise_' : 'self._piecewise_',
        'not' : 'operator.not_', 'not_' : 'operator.not_'
    }

    MathmlToNumpy_symb = {
        'notanumber' : 'numpy.NaN', 'pi' : 'numpy.pi',
        'infinity' : 'numpy.Infinity', 'exponentiale' : 'numpy.e',
        'true' : 'True', 'false' : 'False', 'True' : 'True', 'False' : 'False'
    }

    SymbolReplacements = None
    FunctionReplacements = None

    MathmlToInfix = {
        'and' : 'and', 'or' : 'or', 'true' : 'True', 'false' : 'False', 'xor' : 'xor'
        }

    precedence = (
        ('left','ANDOR'),
        ('left','EQUISYMB'),
        ('left',  'PLUS', 'MINUS'),
        ('left',  'TIMES', 'DIVIDE'),
        ('left',  'POWER'),
        ('right', 'UMINUS')
        )

    # List of token names
    tokens = ('REAL',
              'INT',
              'PLUS',
              'MINUS',
              'TIMES',
              'DIVIDE',
              'POWER',
              'LPAREN',
              'RPAREN',
              'NOTEQUALS',
              'NAME',
              'ANDOR',
              'COMMA',
              'EQUISYMB',
              'PIECEWISE',
              'DELAY')

    def __init__(self):
        self.LexErrors = []

        self.Int = r'\d+'                                      # Integer
        self.Dec = self.Int + '\.' + self.Int                            # Decimal
        self.Exp = r'([E|e][\+|\-]?)' + self.Int                    # Exponent
        self.Real = self.Dec  + '(' + self.Exp + ')?' + '|' + self.Int + self.Exp  # Real - dec w/o optional exp or int with exp

        # Simple tokens
        self.t_REAL = self.Real
        self.t_INT = self.Int
        self.t_PLUS = r'\+'
        self.t_MINUS = r'-'
        self.t_TIMES = r'\*'
        self.t_DIVIDE = r'/'
        self.t_POWER = '\*\*'
        self.t_LPAREN = r'\('
        self.t_RPAREN = r'\)'
        self.t_COMMA = r','
        self.t_NOTEQUALS = r'!='

    def t_NAME(self,t):
        r'numpy\.[\w]*|math\.[\w]*|operator\.[\w]*|self\.[\w]*|[a-zA-Z_][\w]*'

        # names are defined as anything starting with a letter OR numpy. math. or operator.
        t.type = 'NAME'
        # allow self. to be in names but always remove! dodgy testing stage
        t.value = t.value.replace('self.','')
        if t.value == 'and':
            t.type = 'ANDOR'
            t.value = ' %s ' % t.value
        elif t.value == 'or':
            t.type = 'ANDOR'
            t.value = ' %s ' % t.value
        elif t.value == 'xor':
            t.type = 'ANDOR'
            t.value = ' %s ' % t.value
        elif t.value == 'piecewise':
            t.type = 'PIECEWISE'
            t.value = ' %s ' % t.value
        elif t.value == 'delay':
            t.type = 'DELAY'
            t.value = ' %s ' % t.value
        return t

    def t_EQUISYMB(self,t):
        r'>=|<=|!=|==|>|<'
        t.type = 'EQUISYMB'
        t.value = ' %s ' % t.value
        ##  'EQUISYMB', t.value
        return t

    # Define a rule so we can track line numbers
    def t_newline(self,t):
        r'\n+'
        t.lexer.lineno += len(t.value)

    # A string containing ignored characters (spaces and tabs)
    t_ignore  = ' \t'

    # Error handling rule
    def t_error(self,t):
        print("Illegal character '{0}'".format(t.value[0]))
        self.LexErrors.append(t.value[0])
        self.LexOK = False
        t.lexer.skip(1)

    # Build the lexer
    def buildlexer(self,**kwargs):
        # try and find a temporary workspace
        if 'TMP' in os.environ:
            tempDir = os.environ['TMP']
        elif 'TEMP' in os.environ:
            tempDir = os.environ['TEMP']
        else:
            tempDir = os.getcwd()
        os.chdir(tempDir)
        self.lexer = lex.lex(object=self, **kwargs)

    # Test it output
    def testlexer(self,data):
        self.lexer.input(data)
        while 1:
             tok = self.lexer.token()
             if not tok: break
             print(tok)

class MyInfixParser(MyInfixLexer):
    ParseOK = True
    SymbolErrors = None
    ModelUsesNumpyFuncs = 0
    names = None
    functions = None
    output = None
    input = None
    name_prefix = '<pre>'
    name_suffix = '<suf>'
    _runCount = 0
    _runCountmax = 20
    __pwcntr__ = 0
    piecewises = None
    DelayRemoved = False

    def __init__(self):
        MyInfixLexer.__init__(self)
        self.ParseErrors = []
        self.names = []
        self.functions = []
        self.SymbolErrors = []
        self.piecewises = {}

    def setNameStr(self, prefix, suffix):
        self.name_prefix = str(prefix)
        self.name_suffix = str(suffix)

    def p_error(self,t):
        try:
            self.ParseErrors.append(t)
        except:
            print('p_error generated a parsing error')
        tok = yacc.token()
        return tok

    def p_infix(self,t):
        '''Expression : Expression PLUS Expression
                      | Expression MINUS Expression
                      | Expression TIMES Expression
                      | Expression DIVIDE Expression
                      | Expression EQUISYMB Expression
                      | Expression ANDOR Expression
                      | Power
                      | Number
                      | Func
                      | Equivalence
                      | Piecewise
                      | NotEquals
                      | Delay'''
                    # |UMINUS : add if the
                    #  alternative for p_uminus is used

        if len(t.slice)==4:
            t[0] = t[1] + t[2] + t[3]
        else:
            t[0] = t[1]

    def p_notequals(self, t):
        '''NotEquals : NOTEQUALS'''
        t[0] = t[1]

    def p_power(self,t):
        '''Power : Expression POWER Expression'''

        ##  t[0] = 'numpy.power('+ t[1] + ',' + t[3] + ')' #changed to make it DeriVar compatible
        t[0] = 'pow('+ t[1] + ',' + t[3] + ')'

    def p_number(self, t):
        '''Number : REAL
                  | INT
                  | NAME'''
        try:
            t[0] = str(float(t[1]))
        except ValueError:
            if t[1].strip() in self.MathmlToNumpy_symb:
                if self.MathmlToNumpy_symb[t[1]] == None:
                    self.SymbolErrors.append(t[1])
                    print('\nSymbol \"{0}\" not yet supported by PySCeS.'.format(t[1]) )
                    t[0] = 'unknown_symbol_' + t[1]
                else:
                    t[0] = self.MathmlToNumpy_symb[t[1]]
                self.ModelUsesNumpyFuncs = 1
            elif t[1].replace('numpy.','').replace('math.','').replace('operator.','') in self.MathmlToNumpy_symb:
                t[0] = t[1]
            else:
                if self.SymbolReplacements != None and t[1].strip() in self.SymbolReplacements:
                    # replace symb --> prefix.replacement.suffix
                    if self.SymbolReplacements[t[1]] not in self.names:
                        self.names.append(self.SymbolReplacements[t[1]])
                    t[0] = self.name_prefix + self.SymbolReplacements[t[1]] + self.name_suffix
                elif self.FunctionReplacements != None and t[1].strip() in self.FunctionReplacements:
                    # replace symb --> (replacement)
                    t[0] = '(%s)' % self.FunctionReplacements[t[1]]
                else:
                    if t[1] not in self.names:
                        self.names.append(t[1])
                    t[0] = self.name_prefix + t[1] + self.name_suffix

    def p_uminus(self,t):
        '''Expression : MINUS Expression %prec UMINUS'''
        # Alternative '''UMINUS : MINUS Expression'''

        t[0] = t[1] + t[2]

    def p_equivalence(self,t):
        '''Equivalence : ANDOR LPAREN Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                        | ANDOR LPAREN Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression COMMA Expression RPAREN
                       '''

        # this is an almighty hack but i cant see any other way to do it right now ... ALL SUGGESTIONS WELCOME
        # changes and(a,b, .....) or(a,b, .....) to (a and b and ...) (a or b or ...)
        # and not(b) into self._not_(b)
        ##  print 'equivalence', len(t), t[:]
        t[1] = t[1].strip()
        if self.MathmlToInfix.has_key(t[1]):
            t[0] = t[2]
            for tt in range(3, len(t)):
                if t[tt] == ',':
                    if t[1] != 'xor':
                        t[0] += ' %s ' % t[1]
                    else:
                         t[0] += ' %s ' % '!='
                else:
                    t[0] += '%s' % t[tt]

    def p_piecewise(self,t):
        '''Piecewise : PIECEWISE LPAREN ArgListSemiCol RPAREN'''

        t[1] = t[1].strip()
        t[0] = self.MathmlToNumpy_funcs[t[1]] + t[2] + t[3] + t[4]
        pw = t[3].split(';')
        for p in range(len(pw)):
            pw[p] = pw[p].strip()
        name = '__pw%s__' % self.__pwcntr__
        if len(pw) == 3:
            self.__pwcntr__ += 1
            self.piecewises.update({name : {
                                          0 : [pw[1], pw[0]],
                                          'other' : pw[2]
                                          }})
        else:
            self.__pwcntr__ += 1
            self.piecewises.update({name : {}})
            if math.modf(len(pw)/2.0)[0] != 0.0:
                self.piecewises[name].update({'other' : pw.pop(-1)})
            else:
                self.piecewises[name].update({'other' : None})
            for p in range(0,len(pw),2):
                self.piecewises[name].update({p : [pw[p+1], pw[p]]})
        t[0] = self.name_prefix + name + self.name_suffix

    def p_delay(self,t):
        '''Delay : DELAY LPAREN Expression COMMA Expression RPAREN'''

        # for now we just remove the delay on the expression
        self.DelayRemoved = True
        t[0] = t[3]

    def p_function(self,t):
        '''Func : LPAREN ArgList RPAREN
                | NAME LPAREN ArgList RPAREN
                | NAME LPAREN RPAREN
                '''
        # this is to match NAME() which as far as I know is unique to object __calls__
        # as well as differentiate between bracketed functions and expressions:
        # func( S1 ) and ( S/S05 )

        if len(t) == 4:
            if t[1] == '(':
                t[0] = t[1] + t[2] + t[3]
            else:
                t[0] = self.name_prefix + t[1] + t[2] + t[3]
        # convert root(degree,<expr>) to pow(<expr>, 1/degree)
        elif t[1].strip() == 'root':
            t[1] = self.MathmlToNumpy_funcs[t[1]]
            t[3] = '%s, %s' % (t[3][t[3].index(',')+1:], 1.0/float(t[3][:t[3].index(',')])  )
            t[0] = t[1] + t[2] + t[3] + t[4]
        elif t[1].strip() in self.MathmlToNumpy_funcs:
            if self.MathmlToNumpy_funcs[t[1]] == None:
                self.SymbolErrors.append(t[1])
                print('\nFunction \"%s\" not supported by PySCeS' % t[1])
                t[0] = 'unknown_function_'+t[1] + t[2] + t[3] + t[4]
            else:
                try:
                    t[0] = self.MathmlToNumpy_funcs[t[1]] + t[2] + t[3] + t[4]
                except Exception as EX:
                    print('Function Parse error 1 (please report!)\n', EX)
            self.ModelUsesNumpyFuncs = True
        else:
            # t[0] = t[1] + t[2] + t[3]
            # or a numpy fucntion
            if t[1][:6] == 'numpy.' or t[1][:5] == 'math.' or t[1][:9] == 'operator.':
                t[0] = t[1] + t[2] + t[3] + t[4]
            else:
                # assume some arbitrary function definition
                t[0] = self.name_prefix + t[1] + t[2] + t[3] + t[4]
                # add to list of functions
                if t[1] not in self.functions: self.functions.append(t[1])

    # adapted from Andrew Dalke's GardenSnake
    # http://www.dalkescientific.com/writings/diary/GardenSnake.py
    # function arguments f(x,y,z)
    def p_arglist(self,t):
        '''ArgList : Expression
                   | ArgList COMMA Expression'''
        try:
            if len(t) == 2:
                t[0] = t[1]
            elif len(t) == 4:
                t[0] = t[1] + t[2] + t[3]
        except Exception as EX:
            print('Function ArgList error (please report!)\n', EX)

    # expression list f(g(x,y); g(a,b))
    def p_arglist_semicol(self,t):
        '''ArgListSemiCol : Expression
                   | ArgListSemiCol COMMA Expression'''
        try:
            if len(t) == 2:
                t[0] = t[1]
            elif len(t) == 4:
                t[0] = t[1] + '; ' + t[3]
        except Exception as EX:
            print('Function ArgList error (please report!)\n', EX)

    def buildparser(self, **kwargs):
        self.parser = yacc.yacc(module=self, **kwargs)

    def parse(self, data):
        self.ParseErrors = []
        self.LexErrors = []
        self.SymbolErrors = []
        self.names = []
        self.functions = []
        self.input = data
        self.ParseOK = True
        self.LexOK = True
        self.piecewises = {}
        self.DelayRemoved = False
        self.output = self.parser.parse(data)
        ##  assert len(self.SymbolErrors) == 0, '\nUndefined symbols:\n%s' % self.SymbolErrors
        ##  if len(self.SymbolErrors) != 0:
            ##  print '\nUndefined symbols:\n%s' % self.SymbolErrors
        assert self.LexOK, '\nLexer Failure:\n%s' % self.LexErrors
        assert self.ParseOK, '\nParser Failure:\n%s' % self.ParseErrors
        self._runCount += 1
        self.SymbolReplacements = None
        self.FunctionReplacements = None
        if self._runCount > self._runCountmax:
            self._runCount == 0
            # we're back !!!
            imp.reload(lex)
            imp.reload(yacc)
