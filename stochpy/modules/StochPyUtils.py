 #! /usr/bin/env python
"""
StochPy Utils
=============

Module that contains functions that are created by the users of stochpy. New functions will be added in the next releases of stochpy.

Written by T.R. Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net
Last Change: August 06, 2015
"""

from __future__ import division, print_function, absolute_import
import stochpy as _stochpy_
import copy,numpy as np
mod = None

def GetAnalyticalPDF(kon,koff,kdeg,ksyn):
    """ Get the analytical probability density function. The analytical solution is taken from Sharezaei and Swain 2008 - Analytical distributions for stochastic gene expression """    
    import mpmath        
    mpmath.mp.pretty = True
    x_values = np.linspace(0,50,10000)
    y_values = []
    for m in x_values:
        a = ((ksyn/kdeg)**m)*np.exp(-ksyn/kdeg)/mpmath.factorial(m)
        b = mpmath.mp.gamma((kon/kdeg)+m) * mpmath.mp.gamma(kon/kdeg + koff/kdeg)/ (mpmath.mp.gamma(kon/kdeg + koff/kdeg + m)* mpmath.mp.gamma(kon/kdeg))
        c = mpmath.mp.hyp1f1(koff/kdeg,kon/kdeg + koff/kdeg + m,ksyn/kdeg)
        y_values.append(a*b*c)
    return x_values,y_values


def GetAnalyticalWaitingtimes(kon,koff,ksyn):
    """ Get analytical waiting times """
    import mpmath        
    mpmath.mp.pretty = True    
    A = mpmath.sqrt(-4*ksyn*kon+(koff + kon + ksyn)**2)
    x = []
    for i in np.linspace(-20,5,5000):
        x.append(mpmath.exp(i))     
    y = []
    for t in x:
        B = koff + ksyn - (mpmath.exp(t*A)*(koff+ksyn-kon))-kon+A+ mpmath.exp(t*A)*A
        p01diff = mpmath.exp(-0.5*t*(koff + kon + ksyn+A))*B/(2.0*A)
        y.append(p01diff*ksyn)        
    return (x,y) 
    
  
def doSequentialSim(smod,n_generations,cell_division_times):
    """ Model protein synthesis subject to cell division """
    ### Start sequential modelling part ###
    for i in range(1,n_generations):     
        ### divide each species between two daughter cells ###
        for j in range(0,len(smod.data_stochsim.species_labels)): 
            species_amount = smod.SSA.sim_output[-1][1:][j]                 
            if species_amount:
                smod.settings.X_matrix[j] = np.random.binomial(n=species_amount,p=0.5,size=1)
    
        ### replace last time point with species amounts after division ###
        species_after_division = copy.deepcopy(list(smod.settings.X_matrix))
        species_after_division.insert(0,cell_division_times[0:i].sum()) # add time of cell division
        species_after_division.append(np.nan) # no specific reaction occured at cell division
        smod.SSA.sim_output[-1] = copy.deepcopy(species_after_division) 
 
        ### Set settings for new simulation and simulate the next generation ### 
        smod.settings.starttime = copy.deepcopy(smod.SSA.sim_output[-1][0])
        smod.settings.endtime += cell_division_times[i]
        smod.SSA.Execute(smod.settings)
        ### End sequential modelling part ###
    smod.FillDataStochsim()  # add all data to data_stochsim object    
    return smod  

class Utils():
    def __init__(self):
        print("See http://stochpy.sf.net/examples.html for more explanation about these examples")
        global mod
        mod = _stochpy_.SSA()        

    def DoExample1(self):
        """ Immigration Death example (available at http://stochpy.sourceforge.net/examples.html) """        
        mod.Model('ImmigrationDeath.psc')  # Ksyn = 10, Kdeg = 0.2, and mRNA(init) = 50
        lambda_ = 50
        N = 1000000
        data = np.random.poisson(lambda_,N)
        mod.DoStochSim(end=N,mode='steps')
        mod.PlotSpeciesDistributions(linestyle= 'solid')
        n, bins, patches = _stochpy_.plt.hist(data-0.5, max(data)-min(data),normed=1, facecolor='green')
        mod.PrintSpeciesMeans()
        mod.PrintSpeciesStandardDeviations()        
    
    def DoExample2(self):
        """ SBML events and Interpolation example (available at http://stochpy.sourceforge.net/examples.html) """
        mod.Model('dsmts-003-04.xml.psc')
        mod.DoStochSim(end = 50,mode = 'time',trajectories = 1000)
        mod.GetRegularGrid()
        mod.PlotAverageSpeciesTimeSeries()
        mod.PrintAverageSpeciesTimeSeries()
        
    def DoExample3(self):
        """ Burstmodel example (available at http://stochpy.sourceforge.net/examples.html) """
        mod.Model('Burstmodel.psc')  # Parameter values in Burstmodel.psc: kon = koff = 0.05
        ntimesteps = 1000000        
        mod.ChangeParameter("kon",0.05)
        mod.ChangeParameter("koff",0.05)        
        mod.DoStochSim(end=ntimesteps,mode='steps',trajectories=1)
        mod.plot.plotnum = 1 
        mod.PlotSpeciesDistributions(species2plot = 'mRNA', colors=['#00FF00'],linestyle = 'solid') # Usage of html color codes
        mod.PlotWaitingtimesDistributions('R3', colors=['#00FF00'],linestyle = 'None',marker='o')
        
        mod.ChangeParameter("kon",5.0)
        mod.ChangeParameter("koff",5.0)
        mod.DoStochSim(end=ntimesteps,mode='steps',trajectories=1)
        mod.plot.plotnum = 1
        mod.PlotSpeciesDistributions(species2plot = 'mRNA', colors='r',linestyle = 'solid')
        
        kon = 0.05
        koff = 0.05
        kdeg = 2.5
        ksyn = 80.0
        x,y = GetAnalyticalPDF(kon,koff,kdeg,ksyn)
        _stochpy_.plt.figure(1)
        _stochpy_.plt.step(x,y,color ='k')

        kon = 5.0
        koff = 5.0
        x,y = GetAnalyticalPDF(kon,koff,kdeg,ksyn)
        _stochpy_.plt.step(x,y,color ='k')
        _stochpy_.plt.xlabel('mRNA copy number per cell')
        _stochpy_.plt.ylabel('Probability mass')
        _stochpy_.plt.legend(['Bimodal','Unimodal', 'Analytical solution'],numpoints=1,frameon=False)
        _stochpy_.plt.title('')
        _stochpy_.plt.ylim([0,0.045])
        
        mod.PlotWaitingtimesDistributions('R3', colors=['r'],linestyle = 'None',marker='v')
        kon = 0.05
        koff = 0.05
        (x,y) = GetAnalyticalWaitingtimes(kon,koff,ksyn)
        _stochpy_.plt.figure(2)
        _stochpy_.plt.plot(x,y,color ='k')

        kon = 5.0
        koff = 5.0
        (x,y) = GetAnalyticalWaitingtimes(kon,koff,ksyn)
        _stochpy_.plt.plot(x,y,color ='k')
        _stochpy_.plt.xlabel('Time between RNA synthesis events')
        _stochpy_.plt.ylabel('Probability density')
        _stochpy_.plt.legend(['Bimodal','Unimodal', 'Analytical solution'],numpoints=1,frameon=False,loc='lower left')
        _stochpy_.plt.title('')
        _stochpy_.plt.xlim([10**-7,10**3])
        _stochpy_.plt.ylim([10**-9,10**3])
        
    def DoExample4(self):
        """ Second Burstmodel example (available at http://stochpy.sourceforge.net/examples.html) """
        import matplotlib.gridspec as gridspec
        sim_end = 100                
        mod.Model('Burstmodel.psc')
        mod.ChangeParameter("kon",0.05)
        mod.ChangeParameter("koff",0.05)        
        mod.DoStochSim(end=sim_end,mode='time',trajectories=1)
        
        # Use a nice grid to plot 4 figures
        gs = gridspec.GridSpec(4,1,width_ratios=[1],height_ratios=[0.3,1,0.3,1])
        ax1 = _stochpy_.plt.subplot(gs[0])        
        mod.plot.ResetPlotnum()
        mod.PlotSpeciesTimeSeries(species2plot = 'ONstate',IsLegend=False)
        _stochpy_.plt.ion()
        _stochpy_.plt.xlabel('')               # remove xlabel
        _stochpy_.plt.ylabel('')               # remove ylabel       
        _stochpy_.plt.xlim([0,sim_end])        # set x lim
        _stochpy_.plt.xticks([])               # remove x ticks
        _stochpy_.plt.ylim([0,1.1])            # set y lim
        _stochpy_.plt.yticks([])               # remove y lim
        _stochpy_.plt.text(-5.5,0.9,'ON')      
        _stochpy_.plt.text(-5.5,0,'OFF')
        _stochpy_.plt.text(101,0.35,'A',fontsize = 14)

        ax2 = _stochpy_.plt.subplot(gs[1])
        mod.plot.ResetPlotnum()
        mod.PlotSpeciesTimeSeries(species2plot ='mRNA',colors = ['#32CD32'],IsLegend=False,title='')
        _stochpy_.plt.xlim([0,sim_end])
        _stochpy_.plt.xticks([])
        _stochpy_.plt.xlabel('')
        _stochpy_.plt.ylabel('mRNA')
        _stochpy_.plt.yticks([0,20,40,60])
        _stochpy_.plt.text(101,27,'B',fontsize = 14)
        
        mod.ChangeParameter("kon",5.0)
        mod.ChangeParameter("koff",5.0)

        ax3 = _stochpy_.plt.subplot(gs[2])
        mod.plot.ResetPlotnum()
        mod.DoStochSim(end=sim_end,mode='time')
        mod.PlotSpeciesTimeSeries(species2plot ='ONstate',IsLegend=False,title='')
        _stochpy_.plt.xlabel('')
        _stochpy_.plt.ylabel('')
        _stochpy_.plt.xlim([0,sim_end])
        _stochpy_.plt.xticks([])
        _stochpy_.plt.ylim([0,1.1])
        _stochpy_.plt.yticks([])
        _stochpy_.plt.text(-5.5,0.9,'ON')
        _stochpy_.plt.text(-5.5,0,'OFF')
        _stochpy_.plt.text(101,0.35,'C',fontsize = 14)

        ax4 = _stochpy_.plt.subplot(gs[3])
        mod.plot.ResetPlotnum()
        mod.PlotSpeciesTimeSeries(species2plot ='mRNA',colors = ['r'],IsLegend=False,title='')
        _stochpy_.plt.xlim([0,sim_end])        
        _stochpy_.plt.xlabel('Time (min)')
        _stochpy_.plt.ylabel('mRNA')
        _stochpy_.plt.xticks([0,20,40,60,80,100])
        _stochpy_.plt.yticks([0,20,40,60])
        _stochpy_.plt.text(101,27,'D',fontsize = 14)  
        
    def DoExample5(self):
        """ Protein turnover example with and without cell division (available at http://stochpy.sourceforge.net/examples.html) """
        smod = _stochpy_.SSA(model_file='CellDivision.psc',dir=_stochpy_.model_dir) # start SSA module with CellDivision model
        ### do protein synthesis without cell division ###
        smod.DoStochSim(end=10,mode='time')
        smod.PlotSpeciesTimeSeries(species2plot= 'Protein')         
        smod.DoStochSim(end=300,mode='time')
        smod.PlotSpeciesDistributions(species2plot='Protein')

        ### do protein synthesis with cell division ###
        n_generations=10
        cell_division_times = abs(np.random.gamma(scale=0.1,shape=3,size=n_generations))  # cell division times with a gamma distribution
        smod.DoStochSim(end=cell_division_times[0],mode='time',trajectories=1) 
        smod = doSequentialSim(smod,n_generations,cell_division_times)
        smod.PlotSpeciesTimeSeries(species2plot= 'Protein')
        n_generations=500
        cell_division_times = abs(np.random.gamma(scale=0.1,shape=3,size=n_generations))  # cell division times with a gamma distribution
        smod.DoStochSim(end=cell_division_times[0],mode='time',trajectories=1) 
        smod = doSequentialSim(smod,n_generations,cell_division_times)
        smod.PlotSpeciesDistributions(species2plot='Protein')
