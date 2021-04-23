"""
StochPy Version
===============

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: September 25, 2014
"""

MAJOR = 2
MINOR = 3
MICRO = 0
STATUS = ''

def current_version():
    return '%s.%s.%s%s' % (MAJOR, MINOR, MICRO, STATUS)

def current_version_tuple():
    return (MAJOR, MINOR, MICRO)

__version__ = current_version() 
