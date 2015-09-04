"""
Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: Augustus 06, 2014
"""
from __future__ import division, print_function, absolute_import

from . import PyscesInterfaces

class SBML2PSC():
  """
  Module that converts SBML models into PSC models if libxml and libsbml are installed

  Usage:
  >>> converter = stochpy.SBML2PSC()
  >>> converter.SBML2PSC('Burstmodel.xml','/home/user/Stochpy/pscmodels/') 
  """   
  def __init__(self,sbmlfile, sbmldir=None, pscfile=None, pscdir=None,quiet=False):
      """
      Converts a SBML file to a PySCeS MDL input file.

      Input:
       - *sbmlfile*: the SBML file name
       - *sbmldir*: [default = None] the directory of SBML files (if None current working directory is assumed)
       - *pscfile*: [default = None] the output PSC file name (if None *sbmlfile*.psc is used)
       - *pscdir*: [default = None] the PSC output directory (if None the pysces.model_dir is used)
      """      
      PyscesInterfaces.Core2interfaces().convertSBML2PSC(sbmlfile,sbmldir,pscfile,pscdir, netStoich=True,quiet=quiet)
