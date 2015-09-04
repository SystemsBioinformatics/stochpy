"""
PySCeS - Python Simulator for Cellular Systems (http://pysces.sourceforge.net)

Copyright (C) 2004-2014 B.G. Olivier, J.M. Rohwer, J.-H.S Hofmeyr all rights reserved,

Brett G. Olivier (bgoli@users.sourceforge.net)
Triple-J Group for Molecular Cell Physiology
Stellenbosch University, South Africa.

Permission to use, modify, and distribute this software is given under the
terms of the PySceS (BSD style) license. See LICENSE.txt that came with
this distribution for specifics.

NO WARRANTY IS EXPRESSED OR IMPLIED.  USE AT YOUR OWN RISK.
Brett G. Olivier
"""

from __future__ import division, print_function, absolute_import
from .version import __version__

import os
import math, operator
import numpy

from .InfixParser import MyInfixParser

InfixParser = MyInfixParser()
InfixParser.buildlexer()
InfixParser.buildparser(debug=0, debugfile='infix.dbg', tabmodule='infix_tabmodule')
InfixParser.setNameStr('self.', '()')

class MapList(list):
    def __init__(self, *args):
        list.__init__(self,*args)

    def asSet(self):
        return set(self.__getslice__(0, self.__len__()))


class NewCoreBase(object):
    __DEBUG__ = False
    name = None
    annotations = None

    def getName(self):
        return self.name

    def setName(self,name):
        self.name = name

    def get(self, attr):
        """Return an attribute whose name is str(attr)"""
        return self.__getattribute__(attr)

    def getAnnotation(self):
        """Returns an annotation dictionary"""
        if self.annotations == None:
            self.annotations = {}
        return self.annotations.copy()

    def setAnnotation(self, key, value):
        """Set an annotation as a key:value pair"""
        if self.annotations == None:
            self.annotations = {}
        self.annotations.update({key : value})


class NumberBase(NewCoreBase):
    value = None
    value_initial = None

    def __call__(self):
        return self.value

    def getValue(self):
        return self.value

    def setValue(self, v):
        self.value = v


class Compartment(NewCoreBase):
    size = None
    dimensions = None
    Compartment = None
    reactions = None
    species = None
    area = None

    def __init__(self, name, compartment=None):
        self.name = name
        self.Compartment = compartment
        self.reactions = []
        self.species = []

    def __call__(self):
        return self.size

    def setSize(self, size, dim):
        self.size = size
        assert dim in [0,1,2,3], '\nOkeee! {0:d} dimensions?'.format(dim)
        self.dimensions = dim

    def setArea(self, area=None):
        if area == None and self.dimensions == 2:
            self.area = self.size
            if self.__DEBUG__: print('Setting reactive area to size for 2D compartment {0:s}'.format(self.name) )
        elif area == None and self.dimensions == 3:
            self.area = (113.09733552923255*self.size**2.0)**(0.33333333333333331)
            if self.__DEBUG__: print('Setting reactive area to surface area for 3D compartment {0:s} (assuming a sphere geometry)'.format(self.name) )
        self.area = area

    def hasReactions(self):
        return MapList([r.name for r in self.reactions])

    def hasSpecies(self):
        return MapList([s.name for s in self.species])

    def addReaction(self, reaction):
        if reaction.name not in self.hasReactions():
            self.reactions.append(reaction)
            self.__setattr__(reaction.name, reaction)
            if self.__DEBUG__: print('Adding reaction {0:s}'.format(reaction.name) )

    def addSpecies(self, species):
        if species.name not in self.hasSpecies():
            self.species.append(species)
            self.__setattr__(species.name, species)
            if self.__DEBUG__: print('Adding species {0:s}'.format(species.name) )
        else:
            if self.__DEBUG__: print('Species {0:s} already added'.format(species.name) )

    def getDimensions(self):
        return self.dimensions

    def getCompartment(self):
        return self.Compartment

    def hasCompartment(self):
        if self.Compartment != None:
            return True
        else:
            return False

    def isVolume(self):
        if self.dimensions == 3: return True
        else: return False

    def isArea(self):
        if self.dimensions == 2: return True
        else: return False

    def isLength(self):
        if self.dimensions == 1: return True
        else: return False

    def isPoint(self):
        if self.dimensions == 0: return True
        else: return False


class BaseUnit(NewCoreBase):
    '''Base Unit can be of type: time, substance, volume'''
    _types = ('time', 'substance', 'volume','area','length')
    value = 1.0
    type = None

    def __init__(self, name, type):
        self.name = name
        assert type in self._types, '\nType must be one of: {0:s}'.format( str(self._types) )
        self.type = type

    def __call__(self):
        return self.value

    def getType(self):
        return self.type


class SimpleUnit(NewCoreBase):
    exponent = 1.0
    scale = 0.0
    multiplier = 1.0
    baseunit = None
    type = None

    def __init__(self, baseunit, name, exp=1.0, scale=0.0, mult=1.0):
        self.baseunit = baseunit
        self.exponent = exp
        self.scale = scale
        self.multiplier = mult
        self.name = name
        self.type = baseunit.type

    def __call__(self):
        return (self.multiplier*self.baseunit()*10**self.scale)**self.exponent

    def getType(self):
        return self.type


class CompoundUnit(NewCoreBase):
    units = None
    _HAS_USERNAME = False

    def __init__(self, name=None):
        self.units = []
        if name != None:
            self.name = name
            self._HAS_USERNAME = True
        else:
            self.name = ''

    def __call__(self):
        U = 1.0
        for u in self.units:
            U *= u()
        return U

    def addUnit(self, unit):
        self.units.append(unit)
        if not self._HAS_USERNAME:
            self.name = '{0:s}{1}'.format(self.name, unit.getName())

    def getUnits(self):
        return self.units

    def hasUnits(self):
        return MapList([u.getName() for u in self.units])


class Species(NumberBase):
    subs = None
    prods = None
    mods = None
    fixed = False
    Compartment = None
    __amount__ = False

    def __init__(self, name, value):
        self.setName(name)
        self.value = value
        self.value_initial = value
        self.subs = []
        self.prods = []
        self.mods = []

    def getCompartment(self):
        return self.Compartment

    def setCompartment(self, c):
        self.Compartment = c

    def hasCompartment(self):
        if self.Compartment != None:
            return True
        else:
            return False

    def setSubstrate(self, reaction):
        self.__setattr__(reaction.name, reaction)
        self.subs.append(reaction)

    def setProduct(self, reaction):
        self.__setattr__(reaction.name, reaction)
        self.prods.append(reaction)

    def setModifier(self, reaction):
        self.__setattr__(reaction.name, reaction)
        self.mods.append(reaction)

    def isSubstrateOf(self):
        return MapList([r.name for r in self.subs])

    def isProductOf(self):
        return MapList([r.name for r in self.prods])

    def isModifierOf(self):
        return MapList([r.name for r in self.mods])

    def isReagentOf(self):
        return MapList(self.isSubstrateOf() + self.isProductOf())

    def setAmount(self, b):
        self.__amount__ = bool(b)

    def isAmount(self):
        return self.__amount__

class SpeciesAssignmentRule(Species):
    formula = None
    code_string = None
    _names = None
    _functions = None
    type = 'assignment'
    _TIME_ = None

    def __init__(self, name, value):
        Species.__init__(self, name, value)

    def __call__(self):
        exec(self.xcode)
        return self.value

    def addFormula(self, formula):
        formula = formula.replace('self.','')
        self.formula = formula
        InfixParser.setNameStr('self.', '()')
        InfixParser.parse(formula)
        self.code_string = 'self.value={0}'.format(InfixParser.output)
        self._names = InfixParser.names
        self._functions = InfixParser.functions
        self.xcode = compile(self.code_string, '<string>', 'exec')

    def addModelAttr(self, obj):
        self.__setattr__(obj.name, obj)

class Function(NewCoreBase):
    formula = None
    code_string = None
    xcode = None
    value = None
    _names = None
    args = None
    _TIME_ = None

    def __init__(self, name):
        self.setName(name)
        self.args = []

    def __call__(self, *args):
        for ar in range(len(args)):
            self.__setattr__(self.args[ar], args[ar])
        exec(self.xcode)
        return self.value

    def setArg(self, var, value=None):
        self.__setattr__(var, value)
        self.args.append(var)

    def addFormula(self, formula):
        formula = formula.replace('self.','')
        self.formula = formula
        InfixParser.setNameStr('self.', '')
        InfixParser.SymbolReplacements = {'_TIME_':'_TIME_()'}
        InfixParser.parse(formula)
        self._names = InfixParser.names
        self.code_string = 'self.value={0}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, '<string>', 'exec')

class Reaction(NewCoreBase):
    modifiers = None
    substrates = None
    products = None
    stoichiometry = None
    multistoich = None
    multistoich_enabled = False
    parameters = None
    functions = None
    reversible = True
    formula = None
    code_string = None
    rate = None
    xcode = None
    _names = None
    _functions = None
    _TIME_ = None
    Compartment = None

    def __call__(self):
        exec(self.xcode)
        return self.rate

    def __init__(self, name):
        self.setName(name)
        self.modifiers = []
        self.substrates = []
        self.products = []
        self.stoichiometry = {}
        self.parameters = []
        self.functions = []
        self.multistoich = []

    def addSubstrate(self, species):
        self.__setattr__(species.name, species)
        self.substrates.append(species)

    def addProduct(self, species):
        self.__setattr__(species.name, species)
        self.products.append(species)

    def addModifier(self, species):
        self.__setattr__(species.name, species)
        self.modifiers.append(species)

    def addFormula(self, formula):
        formula = formula.replace('self.','')
        self.formula = formula
        InfixParser.setNameStr('self.', '()')
        InfixParser.parse(formula)
        self._names = InfixParser.names
        self._functions = InfixParser.functions
        self.code_string = 'self.rate={0}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, '<string>', 'exec')

    def addParameter(self, par):
        self.__setattr__(par.name, par)
        self.parameters.append(par)

    def addFunction(self, func):
        self.__setattr__(func.name, func)
        self.functions.append(func)

    def hasProducts(self, t=type):
        return MapList([p.name for p in self.products])

    def hasSubstrates(self):
        return MapList([s.name for s in self.substrates])

    def hasModifiers(self):
        return MapList([m.name for m in self.modifiers])

    def hasParameters(self):
        return MapList([p.name for p in self.parameters])

    def hasReagents(self):
        return MapList(self.hasSubstrates() + self.hasProducts())

    def setCompartment(self, compartment):
        self.Compartment = compartment

    def getCompartment(self):
        return self.Compartment

    def hasCompartment(self):
        if self.Compartment != None:
            return True
        else:
            return False


class Parameter(NumberBase):
    association = None

    def __init__(self, name, value):
        self.name = name
        self.value = value
        self.value_initial = value
        self.association = []

    def setAssociation(self, reac):
        self.association.append(reac)
        self.__setattr__(reac.name, reac)

    def isParameterOf(self):
        return MapList([a.name for a in self.association])

class AssignmentRule(Parameter):
    formula = None
    code_string = None
    _names = None
    _functions = None
    type = 'assignment'
    _TIME_ = None
    fixed = False # added so that assignment rules can modify fixed species

    def __init__(self, name, value):
        Parameter.__init__(self, name, value)

    def __call__(self):
        exec(self.xcode)
        return self.value

    def addFormula(self, formula):
        formula = formula.replace('self.','')
        self.formula = formula
        InfixParser.setNameStr('self.', '()')
        InfixParser.parse(formula)
        self.code_string = 'self.value={0}'.format(InfixParser.output)
        self._names = InfixParser.names
        self._functions = InfixParser.functions
        self.xcode = compile(self.code_string, '<string>', 'exec')

    def addModelAttr(self, obj):
        self.__setattr__(obj.name, obj)


class RateRule(NewCoreBase):
    formula = None
    rate = None
    xcode = None
    code_string = None
    _names = None
    _functions = None
    compartment = None

    def __init__(self, name, formula):
        self.name = name
        self.addFormula(formula)

    def __call__(self):
        exec(self.xcode)
        return self.rate

    def addFormula(self, formula):
        formula = formula.replace('self.','')
        self.formula = formula.replace('()','')
        InfixParser.setNameStr('self.', '()')
        InfixParser.parse(self.formula)
        self.code_string = 'self.rate={0}'.format(InfixParser.output)
        self._names = InfixParser.names
        self._functions = InfixParser.functions
        self.xcode = compile(self.code_string, 'RateRule: {0:s}'.format(self.name), 'exec')

    def getFormula(self):
        return self.formula

    def addModelAttr(self, obj):
        self.__setattr__(obj.name, obj)


class ODE(NewCoreBase):
    sdot = None
    value = None
    coefficients = None
    reactions = None
    independent = None
    ode_terms = None
    formula = ''
    formula_alt = ''
    code_string = 'self.value='
    code_string_alt = 'sdot='

    def __init__(self, species, independent=True):
        self.sdot = species
        self.name = 'ODE_'+species.name
        self.reactions = []
        self.coefficients = []
        self.ode_terms = []
        self.independent = independent

    def __call__(self):
        exec(self.code_string)
        return self.value

    def addReaction(self, reaction, coefficient):
        self.reactions.append(reaction)
        self.coefficients.append(coefficient)
        if coefficient > 0.0:
            if coefficient == 1.0:
                term = '+self.{0:s}() '.format(reaction.name)
                aterm = '+({0}) '.format(reaction.code_string.replace('self.rate=',''))
                fterm = '+{0:s} '.format(reaction.name)
                afterm = '+ ({0}) '.format(reaction.formula)
            else:
                term = '+{0}*self.{1}() '.format(abs(coefficient), reaction.name)
                aterm = '+{0}*({1}) '.format(abs(coefficient), reaction.code_string.replace('self.rate=',''))
                fterm = '+{0}*{1}'.format(abs(coefficient), reaction.name)
                afterm = '+ {0}*({1}) '.format(abs(coefficient), reaction.formula)
        else:
            if coefficient == -1.0:
                term = '-self.{0:s}() '.format(reaction.name)
                aterm = '-({0}) '.format(reaction.code_string.replace('self.rate=',''))
                fterm = '-{0:s} '.format(reaction.name)
                afterm = '- ({0}) '.format(reaction.formula)
            else:
                term = '-{0}*self.{1}() '.format(abs(coefficient), reaction.name)
                aterm = '-{0}*({1}) '.format(abs(coefficient), reaction.code_string.replace('self.rate=',''))
                fterm = '-{0}*{1}'.format(abs(coefficient), reaction.name)
                afterm = '- {0}*({1}) '.format(abs(coefficient), reaction.formula)
        self.ode_terms.append(term)
        self.code_string += term
        self.code_string_alt += aterm
        self.formula += fterm
        self.formula_alt += afterm
        self.__setattr__(reaction.name, reaction)

    def hasReactions(self):
        return MapList([r.name for r in self.reactions])

    def getFormula(self):
        return self.code_string

    def getGlobalFormula(self):
        return self.code_string_alt


class StructMatrix(NewCoreBase):
    """
    This class is specifically designed to store structural matrix information
    give it an array and row/col index permutations it can generate its own
    row/col labels given the label src.
    """

    array = None
    ridx = None
    cidx = None
    row = None
    col = None

    def __init__(self, array, ridx, cidx, row=None, col=None):
        """
        Instantiate with array and matching row/col index arrays, optional label arrays
        """
        self.array = array
        self.ridx = ridx
        self.cidx = cidx
        self.row = row
        self.col = col
        self.shape = array.shape

    def __call__(self):
        return self.array

    def getRowsByIdx(self, *args):
        """Return the rows referenced by index (1,3,5)"""
        return self.array.take(args, axis=0)

    def getColsByIdx(self, *args):
        """Return the columns referenced by index (1,3,5)"""
        return self.array.take(args, axis=1)

    def setRow(self, src):
        """
        Assuming that the row index array is a permutation (full/subset)
        of a source label array by supplying that source to setRow it
        maps the row labels to ridx and creates self.row (row label list)
        """
        self.row = [src[r] for r in self.ridx]

    def setCol(self, src):
        """
        Assuming that the col index array is a permutation (full/subset)
        of a source label array by supplying that src to setCol
        maps the row labels to cidx and creates self.col (col label list)
        """
        self.col = [src[c] for c in self.cidx]

    def getRowsByName(self, *args):
        """Return the rows referenced by label ('s','x','d')"""
        assert self.row != None, "\nI need row labels"
        try:
            return self.array.take([self.row.index(l) for l in args], axis=0)
        except Exception as ex:
            print(ex)
            print("\nValid row labels are: {0}".format(self.row) )
            return None

    def getColsByName(self, *args):
        """Return the columns referenced by label ('s','x','d')"""
        assert self.col != None, "\nI need column labels"
        try:
            return self.array.take([self.col.index(l) for l in args], axis=1)
        except Exception as  ex:
            print(ex)
            print("Valid column labels are: {0}".format(self.col) )
            return None

    def getLabels(self, axis='all'):
        """Return the matrix labels ([rows],[cols]) where axis='row'/'col'/'all'"""
        if axis == 'row': return self.row
        elif axis == 'col': return self.col
        else: return self.row, self.col

    def getIndexes(self, axis='all'):
        """Return the matrix indexes ([rows],[cols]) where axis='row'/'col'/'all'"""
        if axis == 'row': return self.ridx
        elif axis == 'col': return self.cidx
        else: return self.ridx, self.cidx

    def getByIdx(self, row, col):
        assert row in self.ridx, '\n{0} is an invalid index'.format(row)
        assert col in self.cidx, '\n{0} is an invalid index'.format(col)
        return self.array[row, col]

    def getByName(self, row, col):
        assert row in self.row, '\n{0} is an invalid name'.format(row)
        assert col in self.col, '\n{0} is an invalid name'.format(col)
        return self.array[self.row.index(row), self.col.index(col)]

    def setByIdx(self, row, col, val):
        assert row in self.ridx, '\n{0} is an invalid index'.format(row)
        assert col in self.cidx, '\n{0} is an invalid index'.format(col)
        self.array[row, col] = val

    def setByName(self, row, col, val):
        assert row in self.row, '\n{0} is an invalid name'.format(row)
        assert col in self.col, '\n{0} is an invalid name'.format(col)
        self.array[self.row.index(row), self.col.index(col)] = val

    def shape(self):
        return self.array.shape

class EventAssignment(NumberBase):
    variable = None
    _names = None
    formula = None
    code_string = None
    xcode = None

    def __call__(self):
        self.variable.value = self.value
        if self.__DEBUG__: print('\tAssigning {0:s} = {1}'.format(self.variable.name, self.value))
        return True

    def __init__(self, name='None'):
        self.setName(name)

    def setVariable(self, var):
        self.variable = var

    def setFormula(self, formula):
        self.formula = formula
        InfixParser.setNameStr('self.', '()')
        ##  InfixParser.SymbolReplacements = {'_TIME_':'_TIME_()'}
        InfixParser.parse(formula)
        self._names = InfixParser.names
        self.code_string = 'self.value={0}'.format(InfixParser.output)
        self.xcode = compile(self.code_string, '<string>', 'exec')
        if self.__DEBUG__: '\t', self.name, self.code_string

    def evaluateAssignment(self):
        exec(self.xcode)


class Event(NewCoreBase):
    trigger = None
    delay = 0.0

    formula = None
    code_string = None
    xcode = None

    state0 = False
    state = False

    assignments = None
    _TIME_ = None
    _ASS_TIME_ = 0.0
    _need_action = False
    _names = None
    _time_symbol = None

    def __init__(self, name):
        self.setName(name)
        self.assignments = []

    def __call__(self, time):
        self._TIME_.set(time)
        exec(self.xcode)
        if self.state0 and not self.state:
            self.state0 = self.state
        if not self.state0 and self.state:
            for ass in self.assignments:
                ass.evaluateAssignment()
            self.state0 = self.state
            self._need_action = True
            self._ASS_TIME_ = self._TIME_() + self.delay
            if self.__DEBUG__: print('event {0:s} is evaluating at {1}'.format(self.name, time))
        if self._need_action and self._TIME_() >= self._ASS_TIME_:
            for ass in self.assignments:
                ass()
            if self.__DEBUG__: print('event {0:s} is assigning at {1} (delay={2})'.format(self.name, time, self.delay))
            self._need_action = False

    def setTrigger(self, formula, delay=0.0):
        self.formula = formula
        self.delay = delay
        InfixParser.setNameStr('self.', '()')
        ##  print self._time_symbol
        if self._time_symbol != None:
            InfixParser.SymbolReplacements = {self._time_symbol : '_TIME_'}
            ##  self.formula = formula.replace(self._time_symbol, '_TIME_')
        InfixParser.parse(formula)
        self._names = InfixParser.names
        self.code_string = 'self.state={0}'.format(InfixParser.output)
        if self._time_symbol != None:
            InfixParser.setNameStr('', '')
            InfixParser.SymbolReplacements = {self._time_symbol : '_TIME_'}
            InfixParser.parse(formula)
            self.formula = InfixParser.output
        self.xcode = compile(self.code_string, '<string>', 'exec')
        if self.__DEBUG__: self.name, self.code_string

    def setTriggerAttributes(self, core):
        # TODO: experimental
        for n in self._names:
            self.__setattr__(n, core.__getattribute__(n))

    def setAssignment(self, var, formula):
        ass = EventAssignment(var.name)
        ass.setVariable(var)
        ass.setFormula(formula)
        self.assignments.append(ass)
        self.__setattr__('_'+var.name, ass)

class PieceWise(NewCoreBase):
    """
    Generic piecewise class written by me!

    - *args* a dictionary of piecewise information generated by the InfixParser
    """
    name = None
    value = None
    formula = None
    code_string = None
    xcode = None
    _names = None
    _TIME_ = None

    def __init__(self, pwd):
        pwd = pwd.copy()
        if pwd['other'] != None:
            other = 'self.value = {0}'.format( pwd.pop('other') )
        else:
            other = 'pass'
            pwd.pop('other')
        InfixParser.setNameStr('self.', '')
        InfixParser.SymbolReplacements = {'_TIME_':'_TIME_()'}
        self._names = []
        if len(list(pwd)) == 1:
            formula = pwd[0][0]
            InfixParser.parse(formula)
            for n in InfixParser.names:
                if n not in self._names and n != '_TIME_()':
                    self._names.append(n)
            formula = InfixParser.output
            self.code_string = 'if {0}:\n    self.value = {1}\nelse:\n    {2}'.format(formula, pwd[0][1], other)
            self.formula = self.code_string.replace('self.','')
        else:
            formula = pwd[0][0]
            InfixParser.parse(formula)
            for n in InfixParser.names:
                if n not in self._names and n != '_TIME_()':
                    self._names.append(n)

            formula = InfixParser.output
            self.code_string = 'if {0}:\n    self.value = {1}\n'.format(formula, pwd[0][1])
            pwd.pop(0)
            for p in pwd:
                formula = pwd[p][0]
                InfixParser.SymbolReplacements = {'_TIME_':'_TIME_()'}
                InfixParser.parse(formula)
                for n in InfixParser.names:
                    if n not in self._names and n != '_TIME_()':
                        self._names.append(n)
                formula = InfixParser.output
                self.code_string += 'elif {0}:\n    self.value = {1}\n'.format(formula, pwd[p][1])
            self.code_string += 'else:\n    {0}'.format(other)
            self.formula = self.code_string.replace('self.','')
        self.xcode = compile(self.code_string, 'PieceWise','exec')


    def __call__(self):
        exec(self.xcode)
        return self.value


class Time(object):
    value = None
    name = '__TIME__'
    def __init__(self, t=0):
        self.value = t

    def __call__(self):
        return self.value

    def set(self, t):
        self.value=t

##  def delay(*args):
    ##  print 'delay() ignored'
    ##  return 1.0

class NewCore(NewCoreBase):
    __nDict__ = None
    reactions = None
    species = None
    species_variable = None
    __model__ = None
    __InitDict__ = None
    __not_inited__ = None
    global_parameters = None
    __parameter_store__ = None
    forcing_functions = None
    __rules__ = None
    __events__ = None
    # new
    __compartments__ = None
    compartments = None
    rate_rules = None
    description = "Pysces Core2"
    __uDict__ = None
    stoichiometric_matrix = None
    struct = None
    ODEs = None
    functions = None
    _TIME_ = None
    events = None
    __sDict__ = None
    __KeyWords__ = None
    __piecewises__ = None
    piecewise_functions = None
    netStoich = None

    def __init__(self, model, iValues=True, netStoich=True):
        # setup core dictionaries
        self.__nDict__ = model.__nDict__
        self.__sDict__ = model.__sDict__
        self.__KeyWords__ = model.__KeyWords__
        if self.__KeyWords__['Modelname'] != None:
            self.setName(self.__KeyWords__['Modelname'])
        else:
            self.setName('PySCeSModel')
        if self.__KeyWords__['Description'] != None:
            self.setDescription(self.__KeyWords__['Description'])
        else:
            self.setDescription('PySCeSModel')

        self.__model__ = model
        self.__InitDict__ = model.__InitDict__
        if not iValues:
            if self.__DEBUG__: print(self.__InitDict__)
            for k in self.__InitDict__.keys():
                self.__InitDict__[k] = getattr(self.__model__, k)
            for c in model.__compartments__:
                model.__compartments__[c]['size'] = getattr(self.__model__, c)
        self.netStoich = netStoich

        self.global_parameters = []
        self.__parameter_store__ = []
        self.__not_inited__ = []
        self.forcing_functions = []
        self.__rules__ = model.__rules__
        self.__uDict__ = model.__uDict__
        self.__piecewises__ = model.__piecewises__
        InfixParser.__pwcntr__ = 0

        # start building objects
        self.__compartments__ = model.__compartments__
        self.addCompartments()
        self._TIME_ = Time()
        self.addPieceWiseFunctions() # this adds any piecewise functions
        self.addSpecies()

        # the order is important from here as eg functions can occur in rate equations
        try:
            self.__functions__ = model.__functions__
        except:
            self.__functions__ = {}
            if self.__DEBUG__: print('No functions')
        self.functions = []
        self.addFunctions()
        self.addReactions()
        self.generateMappings()
        self.setAssignmentRules()
        self.setRateRules()

        # add event support
        self.__events__ = self.__model__.__eDict__
        self.events = []
        self.addEvents()
        self.addPieceWiseFunctions(update=True)  # this updates their attributes

        ##  # get rid of _TIME_ in not intited
        ##  if '_TIME_' in self.__not_inited__:
            ##  self.__not_inited__.pop(self.__not_inited__.index('_TIME_'))
        assert len(self.__not_inited__) < 1, "\nERROR: Uninitialised parameters: {0}".format(self.__not_inited__)

    def __cleanString__(self,s):
        s = s.lstrip()
        s = s.rstrip()
        return s

    def setDescription(self, txt):
        self.description = str(txt)

    def getDescription(self):
        return str(self.description)

    def setGlobalUnits(self, **kwargs):
        for un in kwargs.keys():
            self.__uDict__[un] = (kwargs[un][0], kwargs[un][1])
            if self.__DEBUG__: print("Modified \"{0}\" to be {1}*{0}*10**{2}".format(un, kwargs[un][0], kwargs[un][1]))

    def getGlobalUnits(self):
        return self.__uDict__

    def addPieceWiseFunctions(self, update=False):
        if not update:
            self.piecewise_functions = []
            for pw in self.__piecewises__.keys():
                if self.__DEBUG__: print('Info: adding piecewise function:{0}'.format(pw) )
                P = PieceWise(self.__piecewises__[pw])
                P.setName(pw)
                P.__setattr__('_TIME_', self.__getattribute__('_TIME_'))
                self.piecewise_functions.append(P)
                self.__setattr__(pw, P)
        else:
            for pw in self.piecewise_functions:
                for a in pw._names:
                    pw.__setattr__(a, self.__getattribute__(a))



    def addOneCompartment(self, name, size, dimensions, compartment=None, area=None):
        C = Compartment(name, compartment)
        C.setSize(size, dimensions)
        ##  C.setArea(area)
        self.compartments.append(C)
        self.__setattr__(name, C)

    def addCompartments(self):
        self.compartments = []
        for C in self.__compartments__:
            c2 = self.__compartments__[C]
            if self.__DEBUG__: print('Adding compartment {0}'.format(c2['name']))
            self.addOneCompartment(c2['name'], c2['size'], c2['dimensions'],
                            compartment=c2['compartment'], area=None)

    def addOneSpecies(self, species, value, fix=False, comp=None, amount=False, fullName=None):
        s = Species(species, value)
        ##  if comp != None:
        s.setCompartment(comp)
        s.setAmount(amount)
        s.setAnnotation('sbml_name', fullName)
        if fix: s.fixed = True
        self.__setattr__(species, s)
        self.species.append(s)
        if not fix: self.species_variable.append(s)
        if comp != None:
            comp.addSpecies(s)

    def addSpecies(self):
        self.species = []
        self.species_variable = []
        for s in self.__sDict__:
            name = self.__sDict__[s]['name']
            if s in self.__InitDict__:
                val = self.__InitDict__[s]
            else:
                val = 0.0
            
            fix = self.__sDict__[s]['fixed']
            if self.__sDict__[s]['compartment'] != None:
                comp = self.__getattribute__(self.__sDict__[s]['compartment'])
            else:
                comp = None
            amount = self.__sDict__[s]['isamount']
            fullName = None
            if 'fullName' in self.__sDict__[s]:
                fullName = self.__sDict__[s]['fullName']
            self.addOneSpecies(name, val, fix=fix, comp=comp, amount=amount, fullName=fullName)

    def addOneFunction(self, name, args, formula):
        func = Function(name)
        # TODO: make better
        setattr(func, '_TIME_', self._TIME_)
        for a in args:
            func.setArg(a)
        func.addFormula(formula)
        self.functions.append(func)
        self.__setattr__(name, func)


    def addFunctions(self):
        for f in self.__functions__.keys():
            self.addOneFunction(f,\
                                self.__functions__[f]['args'],\
                                self.__functions__[f]['formula'])

    def addOneReaction(self, rDict):
        r = Reaction(rDict['name'])
        if rDict['compartment'] != None:
            C = self.__getattribute__(rDict['compartment'])
            r.setCompartment(C)
            C.addReaction(r)
        fullName = None
        if 'fullName' in rDict:
            r.setAnnotation('sbml_name', rDict['fullName'])

        # TODO: make better
        setattr(r, '_TIME_', self._TIME_)
        r.addFormula(rDict['RateEq'].replace('self.',''))
        if rDict['Type'] == 'Irrev': r.reversible = False
        # now we can add formulas that occured in the rate equation
        if len(r._functions) > 0:
            for func in r._functions:
                try:
                    r.addFunction(self.__getattribute__(func))
                except Exception as ex:
                    print(ex)
                    print('\nHave you added the function objects yet (addFunctions())')

        #fxnames = self.hasFixedSpecies()
        processed_parameter = []
        # where parameters are defined `locally' per reaction
        for p in rDict['Params']:
            p = p.replace('self.','')
            if p not in self.hasGlobalParameters() and not (p in self.hasFixedSpecies() or p in self.__compartments__):
                if self.__DEBUG__: print("Adding parameter {0} from networkdict".format(p))
                self.addParameter(p)
                par = self.__getattribute__(p)
                par.setAssociation(r)
                r.addParameter(par)
                processed_parameter.append(p)
            elif not (p in self.hasFixedSpecies() or p in self.__compartments__):
                if self.__DEBUG__: print("Updating parameter {0} from networkdict".format(p))
                pidx = self.hasGlobalParameters().index(p)
                self.global_parameters[pidx].setAssociation(r)
                r.addParameter(self.global_parameters[pidx])
                processed_parameter.append(p)

        #print self.hasGlobalParameters()
        # where parameters are not `locally' defined and are extracted from Req (ie from SBML)
        for p in r._names:
            p = p.replace('self.','')
            if p == '_TIME_':
                pass
            elif p in [pw.name for pw in self.piecewise_functions]:
                pass
            elif p in self.hasCompartments() and p not in processed_parameter:
                C = self.__getattribute__(p)
                C.addReaction(r)
                # TODO: this will work until isParameterOf is called on a compartment object
                r.addParameter(C)
                # dirty alternative
                #setattr(r, C.name, C)
                processed_parameter.append(p)
            elif p not in processed_parameter and p not in self.hasGlobalParameters() and p not in self.hasSpecies():
                if self.__DEBUG__: print("Adding parameter {0} from global".format(p) )
                self.addParameter(p)
                par = self.__getattribute__(p)
                par.setAssociation(r)
                r.addParameter(par)
                processed_parameter.append(p)

            elif p not in processed_parameter and p not in self.hasSpecies():
                if self.__DEBUG__: print("Updating parameter {0} from global".format(p) )
                pidx = self.hasGlobalParameters().index(p)
                self.global_parameters[pidx].setAssociation(r)
                r.addParameter(self.global_parameters[pidx])
                processed_parameter.append(p)

        self.__setattr__(rDict['name'], r)
        self.reactions.append(r)

    def addParameter(self, name):
        if not name in self.__piecewises__:
            if name in self.__InitDict__:
                par = Parameter(name, self.__InitDict__[name])
            else:
                par = Parameter(name, 0.0)
                if name not in self.__not_inited__: self.__not_inited__.append(name)
            self.global_parameters.append(par)
            self.__setattr__(name, par)

    def addReactions(self):
        self.reactions = []
        for r in self.__model__.reactions:
            self.addOneReaction(self.__nDict__[r])
        non_parameters = self.hasGlobalParameters()+self.hasSpecies()+self.hasFixedSpecies()
        for k in self.__InitDict__.keys():
            if k not in non_parameters:
                if self.__DEBUG__: print( 'Adding new parameter:', k)
                self.addParameter(k)

    def replaceParameterWithRule(self, ar):
        par = self.__getattribute__(ar.name)
        for r in par.association:
            ar.setAssociation(r)
            setattr(r, ar.name, ar)
            r.parameters[r.hasParameters().index(ar.name)] = ar
        self.global_parameters[self.hasGlobalParameters().index(ar.name)] = ar
        self.__setattr__(ar.name, ar)

    def replaceFixedSpeciesWithRule(self, ar):
        fs = self.__getattribute__(ar.name)
        ar.fixed = fs.fixed
        for r in fs.subs:
            ar.setSubstrate(r)
            setattr(r, ar.name, ar)
            r.substrates[r.hasSubstrates().index(ar.name)] = ar
        for r in fs.prods:
            ar.setProduct(r)
            setattr(r, ar.name, ar)
            r.products[r.hasProducts().index(ar.name)] = ar
        for r in fs.mods:
            ar.setModifier(r)
            setattr(r, ar.name, ar)
            r.modifiers[r.hasModifiers().index(ar.name)] = ar
        self.species[self.hasSpecies().index(ar.name)] = ar
        self.__setattr__(ar.name, ar)

    def replaceSpeciesWithRule(self, ar):
        fs = self.__getattribute__(ar.name)
        for r in fs.subs:
            ar.setSubstrate(r)
            setattr(r, ar.name, ar)
            r.substrates[r.hasSubstrates().index(ar.name)] = ar
        for r in fs.prods:
            ar.setProduct(r)
            setattr(r, ar.name, ar)
            r.products[r.hasProducts().index(ar.name)] = ar
        for r in fs.mods:
            ar.setModifier(r)
            setattr(r, ar.name, ar)
            r.modifiers[r.hasModifiers().index(ar.name)] = ar
        self.species[self.hasSpecies().index(ar.name)] = ar
        self.species_variable[self.hasVariableSpecies().index(ar.name)] = ar
        self.__setattr__(ar.name, ar)

    def setAssignmentRules(self):
        aps = [self.__rules__[ar]['name'] for ar in self.__rules__ if self.__rules__[ar]['type'] == 'assignment']
        ##  for p in self.global_parameters + [self.get(fs) for fs in self.hasFixedSpecies()]:
        for p in self.global_parameters + self.species:
            #print p.name
            if p.name in aps:
                if self.__DEBUG__: print('Assigning: {0:s} = {1}'.format(p.name, self.__rules__[p.name]['formula']))
                p2 = None
                # TODO: make better
                if p.name in self.hasGlobalParameters():
                    p2 = AssignmentRule(p.name, self.__InitDict__[p.name])
                    setattr(p2, '_TIME_', self._TIME_)
                    self.replaceParameterWithRule(p2)
                elif p.name in self.hasFixedSpecies():
                    p2 = SpeciesAssignmentRule(p.name, self.__InitDict__[p.name])
                    p2.setCompartment(p.getCompartment())
                    setattr(p2, '_TIME_', self._TIME_)
                    self.replaceFixedSpeciesWithRule(p2)
                elif p.name in self.hasVariableSpecies():
                    p2 = SpeciesAssignmentRule(p.name, self.__InitDict__[p.name])
                    p2.setCompartment(p.getCompartment())
                    setattr(p2, '_TIME_', self._TIME_)
                    self.replaceSpeciesWithRule(p2)

                assert isinstance(p2, AssignmentRule) or isinstance(p2, SpeciesAssignmentRule), "\nHappy assertion error"
                #print type(p2)
                p2.addFormula(self.__rules__[p.name]['formula'])
                ##  print p2._names
                for n in p2._names+p2._functions:
                    p2.addModelAttr(self.__getattribute__(n))
                ##  # setup initial values
                ##  p2.value_initial = self.p2()
                if p2.name in self.__not_inited__:
                    self.__not_inited__.pop(self.__not_inited__.index(p.name))
        for p in self.global_parameters:
            if p.name in self.hasAssignmentRules():
                # TODO assignment rules need a list of properties
                for ar in p._names:
                    if ar in self.hasAssignmentRules():
                        setattr(p, ar, self.__getattribute__(ar))
        #TODO this is where things will go wrong if fs --> ar contains nested ar's

    def setRateRules(self):
        # TODO mayvbe split into two methods for now read from self.__rules__
        # TODO add functions to rules
        ars = [self.__rules__[ar]['name'] for ar in self.__rules__ if self.__rules__[ar]['type'] == 'rate']
        self.rate_rules = []
        for rr in ars:
            rrobj = RateRule(self.__rules__[rr]['name'], self.__rules__[rr]['formula'])
            ##  print 'RR:', rrobj.name, rrobj._names, rrobj._functions
            for symb in  rrobj._names+rrobj._functions:
                rrobj.addModelAttr(self.__getattribute__(symb))
            self.rate_rules.append(rrobj)
            # TODO investgiate this as it is problematic, the rate rule
            # is not a model property as such more an ODE property
            ##  self.__setattr__(rrobj.name, rrobj)
            if self.__DEBUG__: print( 'Adding RateRule {0} with formula: {1}'.format(rrobj.name, rrobj.formula))

    def addOneEvent(self, e):
        """Add a single event using an event dictionary """
        # translate self.__events__[e] to e

        ev = Event(e['name'])
        ev._time_symbol = e['tsymb']
        ev.setTrigger(e['trigger'], e['delay'])
        # associate model attributes with event
        # TODO: check that this still works
        ev.setTriggerAttributes(self)
        ##  for n in ev._names:
            ##  setattr(ev, n, self.__getattribute__(n))
        # for each assignment
        for ass in e['assignments']:
            ev.setAssignment(self.__getattribute__(ass), e['assignments'][ass])
            assref = getattr(ev, '_'+ass) # don\t like this at all :-(
            # associate model attributes with assignment
            for n in assref._names:
                setattr(assref, n, self.__getattribute__(n))
        self.events.append(ev)
        self.__setattr__(ev.name, ev)
        setattr(ev, '_TIME_', self._TIME_)

    def addEvents(self):
        # TODO: check that you can change the trigger on the fly (might need a setAttr thing in event obj)
        self.events = []
        # for each event
        for e in self.__events__:
            self.addOneEvent(self.__events__[e])

    def generateMappings(self):
        ##  self.netStoich = False
        for reac in self.reactions:
            if self.netStoich:
                for reag in self.__nDict__[reac.name]['Reagents']:
                    if self.__nDict__[reac.name]['Reagents'][reag] < 0.0:
                        reac.addSubstrate(self.__getattribute__(reag.replace('self.','')))
                        self.__getattribute__(reag.replace('self.','')).setSubstrate(self.__getattribute__(reac.name))
                    else:
                        reac.addProduct(self.__getattribute__(reag.replace('self.','')))
                        self.__getattribute__(reag.replace('self.','')).setProduct(self.__getattribute__(reac.name))
                    reac.stoichiometry.setdefault(reag.replace('self.',''), self.__nDict__[reac.name]['Reagents'][reag])
            else:
                for reag in self.__nDict__[reac.name]['AllReagents']:
                    if reag[1] < 0.0:
                        reac.addSubstrate(self.__getattribute__(reag[0].replace('self.','')))
                        self.__getattribute__(reag[0].replace('self.','')).setSubstrate(self.__getattribute__(reac.name))
                    else:
                        reac.addProduct(self.__getattribute__(reag[0].replace('self.','')))
                        self.__getattribute__(reag[0].replace('self.','')).setProduct(self.__getattribute__(reac.name))
                    reac.multistoich.append((reag[0].replace('self.',''), reag[1]))
                    if reag[0].replace('self.','') in reac.stoichiometry:
                        reac.multistoich_enabled = True
                    reac.stoichiometry.setdefault(reag[0].replace('self.',''), reag[1])
            for mod in self.__nDict__[reac.name]['Modifiers']:
                reac.addModifier(self.__getattribute__(mod.replace('self.','')))
                self.__getattribute__(mod.replace('self.','')).setModifier(self.__getattribute__(reac.name))
            ##  print 'I AM LEGEND'
            ##  print reac.stoichiometry
            ##  print reac.multistoich
            ##  print 'reac.multistoich_enabled', reac.multistoich_enabled
            ##  print self.__nDict__[reac.name]['Reagents']
            ##  print self.__nDict__[reac.name]['AllReagents']

    def setStoichiometricMatrix(self):
        vspec = self.hasVariableSpecies()
        react = self.hasReactions()
        nm = numpy.zeros((len(vspec), len(react)),'d')
        for sp in vspec:
            for r in self.get(sp).isReagentOf():
                nm[vspec.index(sp)][react.index(r)] = self.get(r).stoichiometry[sp]
            # this is if absolute stoichiometry value is used
            ##  for r in self.get(sp).isSubstrateOf():
                ##  nm[vspec.index(sp)][react.index(r)] = abs(self.get(r).stoichiometry[sp])
            ##  for r in self.get(sp).isProductOf():
                ##  nm[vspec.index(sp)][react.index(r)] = -abs(self.get(r).stoichiometry[sp])
        self.stoichiometric_matrix = StructMatrix(nm, range(len(vspec)), range(len(react)))
        self.stoichiometric_matrix.setRow(vspec)
        self.stoichiometric_matrix.setCol(react)

    def addODEs(self):
        self.ODEs = []
        for varspec in self.stoichiometric_matrix.row:
            if self.struct != None:
                if varspec not in self.struct.Nr.row:
                    if self.__DEBUG__: print('Creating dependent ODE_{0}'.format(varspec) )
                    ode = ODE(self.get(varspec), independent=False)
                else:
                    if self.__DEBUG__: print('Creating independent ODE_{0}'.format(varspec) )
                    ode = ODE(self.get(varspec), independent=True)
            else:
                if self.__DEBUG__: print( 'Creating independent* ODE_{0} (*assumed - no structural information available)'.format(varspec) )
                ode = ODE(self.get(varspec), independent=True)
            mrow = self.stoichiometric_matrix.getRowsByName(varspec)
            for e in range(len(mrow[0])):
                if mrow[0,e] != 0.0:
                    print('Adding term: {0}*{1}'.format(mrow[0,e], self.stoichiometric_matrix.col[e]))
                    ode.addReaction(self.get(self.stoichiometric_matrix.col[e]), mrow[0,e])
            self.__setattr__(ode.name, ode)
            self.ODEs.append(ode)
            self.__setattr__('xcode_'+ode.name, compile(ode.getGlobalFormula(), '<string>', 'exec'))

    def hasODEs(self):
        return MapList([o.name for o in self.ODEs])

    def evalODEs(self, odes):
        return [v() for v in odes]

    def evalXcode(self, ode):
        exec(self.__getattribute__('xcode_'+ode.name))
        return sdot

    def hasFunctions(self):
        return MapList([f.name for f in self.functions])

    def hasReactions(self):
        return MapList([r.name for r in self.reactions])

    def hasSpecies(self):
        return MapList([s.name for s in self.species])

    def hasFixedSpecies(self):
        return MapList([s.name for s in self.species if s.fixed])

    def hasVariableSpecies(self):
        return MapList([s.name for s in self.species if not s.fixed])

    def findReactionsThatIncludeAllSpecifiedReagents(self, *args):
        assert len(args) > 1, '\nNeed two or more species for this one!'
        setlist = [self.__getattribute__(s).isReagentOf().asSet() for s in args]
        isect = setlist[0]
        for s in setlist:
            isect.intersection_update(s)
        return MapList(isect)

    def hasGlobalParameters(self):
        return MapList(p.name for p in self.global_parameters)

    def hasAssignmentRules(self):
        return MapList([ar.name for ar in self.global_parameters+self.species if hasattr(ar, 'type')=='assignemnt'])

    def hasAssignmentRules(self):
        return MapList([ar.name for ar in self.global_parameters+self.species if hasattr(ar, 'type')=='rate'])

    def hasEvents(self):
        return MapList(e.name for e in self.events)

    def hasCompartments(self):
        return MapList(c.name for c in self.compartments)
