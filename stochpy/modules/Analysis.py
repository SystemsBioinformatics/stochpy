#! /usr/bin/env python
"""
Analysis
========

This module provides functions for Stochastic Simulation Algorithms Analysis (SSA). Implemented SSAs import this module to perform their analysis. Plotting of time series species, propensities), distributions (species, propensities, distributions), autocorrelations, and autocovariances (species, propensities) is possible.

Written by TR Maarleveld, Amsterdam, The Netherlands
E-mail: tmd200@users.sourceforge.net

"""

from __future__ import division, print_function, absolute_import

from stochpy import _IsPlotting
if _IsPlotting:
    from stochpy import plt
    from stochpy import matplotlib
    from matplotlib import gridspec,colors as clr

from stochpy import _IsNumPy
if _IsNumPy:
    import numpy as np
else:
    sys.exit()

import copy,sys

def getDataForTimeSimPlot(Arr_data,npoints = 100000,quiet=False):
    """
    Input:
     - *Arr_data* (numpy array)
     - *npoints* [default = 10000] (integer)
    """
    len_data = len(Arr_data)
    if (len_data > npoints): # use npoints only if datasets become too large
        L_data2plot = [Arr_data[0]]
        step_size = len_data//int(abs(npoints))
        for i in range(step_size,len_data,step_size):
            t = Arr_data[i][0]
            data_point = copy.deepcopy(L_data2plot[-1][1:].tolist())
            data_point.insert(0,t)
            L_data2plot.append(data_point)
            L_data2plot.append(Arr_data[i])
        if not quiet:
            print("Info: Plotting {0:d} out of {1:d} points. Use the argument 'npoints' to alter the number of plotted events.".format(npoints,len_data) )
    else:
        L_data2plot = copy.deepcopy(Arr_data.tolist())
        j=1
        for i in range(1,len_data):
            t = Arr_data[i][0]
            data_prev = copy.deepcopy(Arr_data[i-1]) # data of previous ...
            data_prev[0] = t
            L_data2plot.insert(j,data_prev)
            j+=2
    return np.array(L_data2plot)


def Count(data,edges):
    """
    Input:
     - *data* (list)
     - *edges* (list)
    """
    n_edges = len(edges)
    L_output = np.zeros(n_edges)
    for value in data:
        for i in range(n_edges-1):
            if (value >= edges[i]) and (value < edges[i+1]):
                L_output[i]+=1
    return np.array(L_output)


def GetSpeciesDistributions(sim_output,species):
    """
    Get distributions, means, standard deviations, and the (raw) moments

    Input:
    - *sim_output* (list)
    - *species* (list)

    Mean = mu = sum(x*P(x))
    Variance = sum(x^2 * p(x)) - mu**2

    Output:
    - *L_probability_mass*
    - *D_means*
    - *D_stds*
    - *D_moments*
    """
    n_species = len(species)
    L_distributions = [{} for i in range(n_species)]
    starttime = sim_output[0][0]
    endtime = sim_output[-1][0]
    n_datapoints = len(sim_output)
    D_means = {}
    D_stds = {}
    D_moments = {}
    L_probability_mass = []
    if n_datapoints > 1:
        for t in range(n_datapoints-1):
            for i in range(n_species):
                try:
                    L_distributions[i][int(sim_output[t][i+1])] += sim_output[t+1][0] - sim_output[t][0]
                except KeyError:
                    L_distributions[i][int(sim_output[t][i+1])] = sim_output[t+1][0] - sim_output[t][0]
        for i,s_id in enumerate(species):
            x = np.array(sorted(L_distributions[i]),dtype=int)
            p_x = np.array([L_distributions[i][x_i] for x_i in x])/float(endtime-starttime) # probability = dt/T

            mu = (x*p_x).sum()
            mu_sq = (x**2*p_x).sum()
            var = mu_sq - mu**2
            std = var**0.5
            L_probability_mass.append([x,p_x])

            D_means[s_id] = mu
            D_stds[s_id] = std

            D_moments[s_id] = {}
            D_moments[s_id]['1'] = mu
            D_moments[s_id]['2'] = mu_sq
            D_moments[s_id]['3'] = (x**3*p_x).sum()
            D_moments[s_id]['4'] = (x**4*p_x).sum()

    return (L_probability_mass,D_means,D_stds,D_moments)


def GetDataDistributions(sim_output,identifiers):
    """
    Get distributions, means, standard deviations, and the (raw) moments

    This function is different, because it does not assume integers, like GetSpeciesDistributions()

    Input:
    - *sim_output* (list)
    - *identifiers* (list)

    Mean = mu = sum(x*P(x))
    Variance = sum(x^2 * p(x)) - mu**2

    Output:
    - *L_probability_mass*
    - *D_means*
    - *D_stds*
    - *D_moments*
    """
    n_identifiers = len(identifiers)
    L_distributions = [{} for i in range(n_identifiers)]
    starttime = sim_output[0][0]
    endtime = sim_output[-1][0]
    n_datapoints = len(sim_output)
    D_means = {}
    D_stds = {}
    D_moments = {}
    L_probability_mass = []
    if n_datapoints > 1:
        for t in range(n_datapoints-1):
            for i in range(n_identifiers):
                try:
                    L_distributions[i][sim_output[t][i+1]] += sim_output[t+1][0] - sim_output[t][0]
                except KeyError:
                    L_distributions[i][sim_output[t][i+1]] = sim_output[t+1][0] - sim_output[t][0]
        for i,id in enumerate(identifiers):
            x = np.array(sorted(L_distributions[i]))
            p_x = np.array([L_distributions[i][x_i] for x_i in x])/float(endtime-starttime) # probability = dt/T

            mu = (x*p_x).sum()
            mu_sq = (x**2*p_x).sum()
            var = mu_sq - mu**2
            std = var**0.5
            L_probability_mass.append([x,p_x])

            D_means[id] = mu
            D_stds[id] = std

            D_moments[id] = {}
            D_moments[id]['1'] = mu
            D_moments[id]['2'] = mu_sq
            D_moments[id]['3'] = (x**3*p_x).sum()
            D_moments[id]['4'] = (x**4*p_x).sum()

    return (L_probability_mass,D_means,D_stds,D_moments)


def LogBin(data,factor):
    """
    Function that creates log bins

    Input:
     - *data* (list)
     - *factor* (float) determines the width of the bins
    Output:
     - *L_x* (list)
     - *L_y* (list)
     - *nbins* (integer)
    """
    xmin = float(min(data))
    nbins = int(np.ceil(np.log(max(data)/xmin)/np.log(factor)))
    L_x = None
    L_y = None
    if nbins:
        L_edges = np.zeros(nbins)
        L_edges[0] = xmin
        for i in range(1,nbins): # 1,nbins
            L_edges[i] = L_edges[i-1]*factor

        L_x  = L_edges[0:(nbins-1)]+np.diff(L_edges)/2
        L_dp = Count(data,L_edges)
        L_ry = np.array(L_dp[0:(nbins-1)])
        L_dedges = np.array(np.diff(L_edges))
        L_y = L_ry/(sum(L_ry)*L_dedges)
    return(L_x,L_y,nbins)


def ObtainWaitingtimes(data_stochsim,reactions):
    """
    This function extracts the waiting times for each reaction of the model from the used SSA output.

    Input:
     - *data_stochsim* (python data object) that stores all simulation data
     - *reactions* (list)
    output:
     - *D_waiting_times* (dict)

    Note: It is impossible to use this function in combination with the Tau-leaping method, because the Tau-Leaping results are not exact!
    """
    L_time = data_stochsim.time.flatten()
    L_fired_reactions = data_stochsim.fired_reactions              # Reactions that fired at some time point
    D_waiting_times = {}
    D_last_time_fired = {}
    nreactions = len(reactions)
    for r_id in reactions:
        D_waiting_times[r_id] = []                                 # create a list that will contain event waiting times for reaction r

    for (current_time,r_index) in zip(L_time[1:],L_fired_reactions[1:]): # Updated Oktober 1st
        for i in range(1,nreactions+1):                            # fired reactions are (1,2,3, .... nreactions)
            if r_index == i:
                if r_index in D_last_time_fired:
                    r_name = reactions[int(r_index-1)]
                    D_waiting_times[r_name].append(current_time - D_last_time_fired[r_index]) # Add inter-arrival time
                    D_last_time_fired[r_index] = current_time      # Update last firing time
                else:
                    D_last_time_fired[r_index] = current_time      # Initial firing time

            elif r_index == -i:                                    # Handle delayed completions 01-10-2014
                r_name_compl = reactions[ int(abs(r_index)-1) ] + '_Completion'
                if r_index in D_last_time_fired:
                    D_waiting_times[r_name_compl].append(current_time - D_last_time_fired[r_index]) # Add inter-arrival time
                    D_last_time_fired[r_index] = current_time      # Update last firing time
                else:
                    D_last_time_fired[r_index] = current_time      # Initial firing time
                    D_waiting_times.setdefault(r_name_compl, [])   # Set keyname if not present\
        #print current_time,D_last_time_fired
    return D_waiting_times


def GetAverageResults(regular_grid):
    """
    Gets the averaged output of multiple trajectories

    Input:
     - *regular_grid* (nested list)
    """
    means = []
    stds = []
    for data in regular_grid:
        means.append(np.mean(data,0))
        stds.append(np.std(data,0))
    return (np.array(means).transpose(),np.array(stds).transpose()) # test: 27 july 15


def RemoveBias(x, axis):
    "Subtracts an estimate of the mean from signal x at axis"
    padded_slice = [slice(d) for d in x.shape]
    padded_slice[axis] = np.newaxis
    mn = np.mean(x, axis=axis)
    return x - mn[tuple(padded_slice)]


def AutoCov(s, **kwargs):
    """
    Returns the autocovariance of signal s at all lags.

    Notes:
    Adheres to the definition
    sxx[k] = E{S[n]S[n+k]} = cov{S[n],S[n+k]}
    where E{} is the expectation operator, and S is a zero mean process
    """
    # only remove the mean once, if needed
    debias = kwargs.pop('debias', True)
    axis = kwargs.get('axis', -1)
    if debias:
        s = RemoveBias(s, axis)
    kwargs['debias'] = False
    return CrossCov(s, s, **kwargs)


def FFTconvolve(in1, in2, mode="full", axis=None):
    """ Convolve two N-dimensional arrays using FFT. See convolve. """
    s1 = np.array(in1.shape)
    s2 = np.array(in2.shape)
    complex_result = (np.issubdtype(in1.dtype, np.complex) or
                      np.issubdtype(in2.dtype, np.complex))
    if axis is None:
        size = s1+s2-1
        fslice = tuple([slice(0, int(sz)) for sz in size])
    else:
        equal_shapes = s1==s2
        # allow equal_shapes[axis] to be False
        equal_shapes[axis] = True
        assert equal_shapes.all(), 'Shape mismatch on non-convolving axes'
        size = s1[axis]+s2[axis]-1
        fslice = [slice(l) for l in s1]
        fslice[axis] = slice(0, int(size))
        fslice = tuple(fslice)

    # Always use 2**n-sized FFT
    fsize = int(2**np.ceil(np.log2(size)))
    if axis is None:
        IN1 = np.fft.fftpack.fftn(in1,fsize)
        IN1 *= np.fft.fftpack.fftn(in2,fsize)
        ret = np.fft.fftpack.ifftn(IN1)[fslice].copy()
    else:
        IN1 = np.fft.fftpack.fft(in1,fsize,axis=axis)
        IN1 *= np.fft.fftpack.fft(in2,fsize,axis=axis)
        ret = np.fft.fftpack.ifft(IN1,axis=axis)[fslice].copy()
    del IN1
    if not complex_result:
        ret = ret.real
    if mode == "full":
        return ret
    elif mode == "same":
        if np.product(s1,axis=0) > np.product(s2,axis=0):
            osize = s1
        else:
            osize = s2
        return _centered(ret,osize)
    elif mode == "valid":
        return _centered(ret,abs(s2-s1)+1)


def CrossCov(x, y, axis=-1, all_lags=False, debias=True):
    """
    Returns the crosscovariance sequence between two ndarrays.
    This is performed by calling fftconvolve on x, y[::-1]

    Input:
    - *x*: ndarray
    - *y*: ndarray
    - *axis*: time axis
    - *all_lags*: {True/False}
       whether to return all nonzero lags, or to clip the length of s_xy
       to be the length of x and y. If False, then the zero lag covariance
       is at index 0. Otherwise, it is found at (len(x) + len(y) - 1)/2
    - *debias*: {True/False}
       Always removes an estimate of the mean along the axis, unless
       told not to.

    Notes:
    cross covariance is defined as
    sxy[k] := E{X[t]*Y[t+k]}, where X,Y are zero mean random processes
    """
    if x.shape[axis] != y.shape[axis]:
        raise ValueError('CrossCov() only works on same-length sequences for now')
    if debias:
        x = RemoveBias(x, axis)
        y = RemoveBias(y, axis)
    slicing = [slice(d) for d in x.shape]
    slicing[axis] = slice(None,None,-1)
    sxy = FFTconvolve(x, y[tuple(slicing)], axis=axis, mode='full')
    N = x.shape[axis]
    sxy /= N
    if all_lags:
        return sxy
    slicing[axis] = slice(N-1,2*N-1)
    return sxy[tuple(slicing)]


def Autocorrelation(s, **kwargs):
    """
    Returns the autocorrelation of signal s at all lags.

    Notes:
    Adheres to the definition
    rxx[k] = E{S[n]S[n+k]}/E{S*S} = cov{S[n],S[n+k]}/sigma**2
    where E{} is the expectation operator, and S is a zero mean process
    """
    # only remove the mean once, if needed
    debias = kwargs.pop('debias', True)
    axis = kwargs.get('axis', -1)
    if debias:
        s = RemoveBias(s, axis)
        kwargs['debias'] = False
    sxx = AutoCov(s, **kwargs)
    all_lags = kwargs.get('all_lags', False)
    if all_lags:
        i = (2*s.shape[axis]-1)/2
        sxx_0 = sxx[i]
    else:
        sxx_0 = sxx[0]
    if not sxx_0:
        sxx = [np.nan for i in range(len(sxx))] # Modification
    else:
        sxx /= sxx_0
    return sxx


class DoPlotting():
    """
    This class initiates the plotting options.

    Input:
     - *species_labels* (list) [S1,S2, ..., Sn]
     - *rate_labels* (list) [R1, R2, ..., Rm]
    """
    def __init__(self,species_labels,rate_labels,plotnum=1,quiet = False):
        self.species_labels = species_labels
        self.rate_labels = rate_labels
        self.number_of_rates = len(rate_labels)
        self.plotnum  = plotnum
        # https://github.com/matplotlib/matplotlib/blob/master/lib/matplotlib/colors.py
        self.colors = ['#0000FF','#00CC00','#FF0033','#FF00CC','#6600FF','#FFFF00','#000000','#CCCCCC','#00CCFF','#99CC33','#FF6666', '#FF99CC','#CC6600','#003300','#CCFFFF','#9900FF','#CC6633','#FFD700','#C0C0C0']
        self.quiet = quiet


    def ResetPlotnum(self):
        """ Reset figure numbers if trajectories > 1 """
        self.plotnum = 1


    def TimeSeries(self,data,npoints,datatype,labels,trajectory_index,linestyle,linewidth,marker,colors,title,xlabel,ylabel,is_legend,legend_location):
        """
        Tracks the propensities and/or species over time.

        Input:
         - *data* (array)
         - *npoints* (integer)
         - *datatype* (list)
         - *labels* (list)
         - *trajectory_index* (integer)
         - *linestyle* (string)
         - *linewidth* (float)
         - *title* (string)
         - *xlabel* (string)
         - *ylabel* (string)
         - *is_legend* (boolean)
        """
        plt.figure(self.plotnum)
        datatype_indices = [labels.index(Id) for Id in datatype]

        data = getDataForTimeSimPlot(data,npoints,self.quiet)

        Arr_time = data[:,0]
        if len(datatype) == 1:
            j = trajectory_index
        else:
            j=0

        for i in datatype_indices:
            y = data[:,i+1]
            if colors == None:
                if j >= len(self.colors):
                    j=0
            elif isinstance(colors,list):
                if j >= len(colors):
                    j=0
            elif isinstance(colors,str):
                colors = [colors]
                j=0

            if colors == None:
                if marker == '' and linestyle == 'solid':
                    plt.plot(Arr_time,y, ls = linestyle,lw = linewidth,color = self.colors[j])
                else:
                    plt.plot(Arr_time,y,marker,ls = linestyle,lw = linewidth,color = self.colors[j])
            else:
                if clr.is_color_like(colors[j]):
                    plt.plot(Arr_time,y,marker,ls = linestyle,lw = linewidth,color = colors[j])
                else:
                    print("*** WARNING ***: '{0}' is not recognized as a valid color code".format(colors[j]) )
                    plt.plot(Arr_time,y,marker,ls = linestyle,lw = linewidth,color = self.colors[j])
                    colors = None
            j+=1
        if is_legend:
            plt.legend(datatype,numpoints=1,frameon=True,loc=legend_location)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)


    def Autocorrelations(self,lags,data,datatype,labels,trajectory_index,linestyle,linewidth,marker,colors,title,xlabel,ylabel,is_legend,legend_location):
        """
        Input:
         - *lags*
         - *data* (array)
         - *datatype* (list)
         - *labels* (list)
         - *trajectory_index* (integer)
         - *linestyle* (string)
         - *linewidth* (float)
         - *marker* string)
         - *colors* (list)
         - *title* (string)
         - *xlabel* (string)
         - *ylabel* (string)
         - *is_legend* (boolean)
        """
        plt.figure(self.plotnum)
        datatype_indices = [labels.index(Id) for Id in datatype]
        if len(datatype) == 1:
            j = trajectory_index
        else:
            j=0

        for i in datatype_indices:
            if colors == None:
                if j >= len(self.colors):
                    j=0
            elif isinstance(colors,list):
                if j >= len(colors):
                    j=0
            elif isinstance(colors,str):
                colors = [colors]
                j=0

            y = data[i][0:len(lags)]
            if colors == None:
                plt.plot(lags,y,marker,ls = linestyle,lw = linewidth, color = self.colors[j])
            else:
                if clr.is_color_like(colors[j]):
                    plt.plot(lags,y,marker,ls = linestyle,lw = linewidth, color = colors[j])
                else:
                    print("*** WARNING ***: '{0}' is not recognized as a valid color code".format(colors[j]) )
                    plt.plot(lags,y,marker,ls = linestyle,lw = linewidth, color = self.colors[j])
                    colors = None
            j+=1
        if is_legend:
            plt.legend(datatype,numpoints=1,frameon=True,loc=legend_location)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)


    def Distributions(self,distributions,datatype,labels,trajectory_index,linestyle,linewidth,colors,title,xlabel,ylabel,is_legend=True,legend_location='upper right',bin_size=1,histtype = 'step',orientation='vertical',multiplotting=False):
        """
        Plots the distributions of species and/or propensities

        density=False because the total probability is determined by summation not by integration.

        Input:
         - *distributions* (nested list)
         - *datatype* (list)
         - *labels* (list)
         - *trajectory_index* (integer)
         - *linestyle* (string)
         - *linewidth* (float)
         - *colors* (list)
         - *title* (string)
         - *xlabel* (string)
         - *ylabel* (string)
         - *is_legend* (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
         - *bin_size* (string) [default = 1]
         - *histtype* (string)) [default = 'step']
         - *orientation* (string) [default = 'vertical']
         - *multiplotting* (boolean) [default = False]
        """
        plt.figure(self.plotnum)
        datatype_indices = [labels.index(Id) for Id in datatype]
        if len(datatype) == 1:
            j = trajectory_index
        else:
            j=0

        for i in datatype_indices:
            dat_min = distributions[i][0].min()
            dat_max = distributions[i][0].max()
            n_bins = 1 + (dat_max-dat_min) / bin_size # Just take one trajectory as reference
            L_bin_edges = np.linspace(dat_min-bin_size/2.0,dat_max+bin_size/2.0,int(n_bins+1))

            if colors == None:
                if j >= len(self.colors):
                    j=0
            elif isinstance(colors,list):
                if j >= len(colors):
                    j=0
            elif isinstance(colors,str):
                colors = [colors]
                j=0

            if colors == None:
                print('#'*20)
                output = plt.hist(distributions[i][0], L_bin_edges, weights = distributions[i][1], ls = linestyle, lw = linewidth, color = self.colors[j], histtype = histtype, orientation=orientation, )
                print('just ran this line')
                output = plt.hist(distributions[i][0], L_bin_edges, weights = distributions[i][1], ls = linestyle, lw = linewidth, color = self.colors[j], histtype = histtype, orientation=orientation, density=False)
                print('just ran this line')
                print('#'*20)

            else:
               if clr.is_color_like(colors[j]):
                    output = plt.hist(distributions[i][0],L_bin_edges,weights = distributions[i][1],ls = linestyle,lw = linewidth,color = colors[j],histtype = histtype,orientation=orientation,density=False)
               else:
                    print("*** WARNING ***: '{0}' is not recognized as a valid color code".format(colors[j]) )
                    output = plt.hist(distributions[i][0],L_bin_edges,weights = distributions[i][1],ls = linestyle,lw = linewidth,color = self.colors[j],histtype = histtype,orientation=orientation,density=False)
                    colors = None
            j+=1
        if is_legend:
            plt.legend(datatype,numpoints=1,frameon=True,loc=legend_location)
        plt.title(title)
        if orientation.lower() == 'horizontal':
            plt.xlabel(ylabel)
            plt.ylabel(xlabel)
            if multiplotting:
                plt.xticks([0,max(output[0])*1.1])
        else:
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            if multiplotting:
                plt.yticks([0,max(output[0])*1.1])


    def WaitingtimesDistributions(self,waiting_times,rates,trajectory_index,linestyle,linewidth, marker,colors,title,xlabel,ylabel,is_legend,legend_location):
        """
        Plots the waiting times for each reaction in the model.
        Makes use of ObtainWaitingtimes to derive the waiting times out of the SSA output.

        Input:
         - *waiting_times* (dict)
         - *rates* (list)
         - *trajectory_index* (integer)
         - *linestyle* (string)
         - *linewith* (float)
         - *marker* (string)
         - *colors* (list)
         - *title* (string)
         - *xlabel* (string)
         - *ylabel* (string)
         - *is_legend* (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        plt.figure(self.plotnum)
        if len(rates) == 1:
            j = trajectory_index
        else:
            j=0

        L_legend_names = []
        for r_id in rates:
            L_waiting_times = waiting_times[r_id]         # get list of waiting times for a given reaction
            if len(L_waiting_times) > 1:			      # At least 2 waiting times are necessary per reaction
                (x,y,nbins) = LogBin(L_waiting_times,1.5) # Create logarithmic bins (HARDCODED 1.5)

                if x is not None:
                    if colors == None:
                        if j >= len(self.colors):
                            j=0
                    elif isinstance(colors,list):
                        if j >= len(colors):
                            j=0
                    elif isinstance(colors,str):
                        colors = [colors]
                        j=0
                    if colors == None:
                        plt.loglog(x,y,marker,ls = linestyle,lw=linewidth,color = self.colors[j])
                    else:
                        if clr.is_color_like(colors[j]):
                            plt.loglog(x,y,marker,ls = linestyle,lw=linewidth,color = colors[j])
                        else:
                            print("*** WARNING ***: '{0}' is not recognized as a valid color code".format(colors[j]) )
                            plt.loglog(x,y,marker,ls = linestyle,lw=linewidth,color = self.colors[j])
                            colors = None
                    L_legend_names.append(r_id)
                    j+=1
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        if is_legend:
            plt.legend(L_legend_names,numpoints=1,frameon=True,loc=legend_location)


    def AverageTimeSeries(self,means,stds,time,nstd,datatype,labels,linestyle,linewidth,marker,ms,colors,title,xlabel,ylabel,is_legend,legend_location):
        """
        Plots the average and standard deviation of datatype on a regular grid.

        Input:
         - *means* (array)
         - *stds* (array)
         - *time* (array)
         - *nstd* (float)
         - *datatype* (list)
         - *labels* (list)
         - *linestyle* (string)
         - *linewidth* (float)
         - *marker* (string)
         - *ms* (float)
         - *colors* (list)
         - *title* (string)
         - *xlabel* (string)
         - *ylabel* (string)
         - *is_legend* (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        assert nstd > 0, "Error: The number of STDs must be a value larger than zero"
        plt.figure(self.plotnum)
        datatype_indices = [labels.index(Id) for Id in datatype]
        j=0
        for i in datatype_indices:
            if colors == None:
                if j >= len(self.colors):
                    j=0
            elif isinstance(colors,list):
                if j >= len(colors):
                    j=0
            elif isinstance(colors,str):
                colors = [colors]
                j=0

            # plot with y-axis error bars
            if colors == None:
                plt.errorbar(time,means[:,i],yerr = nstd*np.array(stds[:,i]),color = self.colors[j],ls = linestyle,lw=linewidth,marker = marker,ms=ms,label = labels[i])
            else:
                if clr.is_color_like(colors[j]):
                    plt.errorbar(time,means[:,i],yerr = nstd*np.array(stds[:,i]),color = colors[j],ls = linestyle,lw=linewidth,marker = marker,ms=ms,label = labels[i])
                else:
                    print("*** WARNING ***: '{0}' is not recognized as a valid color code".format(colors[j]) )
                    plt.errorbar(time,means[:,i],yerr = nstd*np.array(stds[:,i]),color = self.colors[j],ls = linestyle,lw=linewidth,marker = marker,ms=ms,label = labels[i])
                    colors = None
            j+=1
        if is_legend:
            plt.legend(numpoints=1,frameon=True,loc=legend_location)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)


    def AverageDistributions(self,means,stds,nstd,datatype,labels,linestyle,linewidth,marker,colors,title,xlabel,ylabel,is_legend,legend_location):
        """
        Plots the average and standard deviation.

        Input:
         - *means* (nested list)
         - *stds* (nested list)
         - *nstd* (float)
         - *labels* (list)
         - *linestyle* (string)
         - *linewidth* (float)
         - *marker* (string)
         - *colors* (list)
         - *title* (string)
         - *xlabel* (string)
         - *ylabel* (string)
         - *is_legend* (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        assert nstd > 0, "Error: The number of STDs must be a value larger than zero"
        plt.figure(self.plotnum)
        datatype_indices = [labels.index(Id) for Id in datatype]
        j=0
        for i in datatype_indices:
            if colors == None:
                if j >= len(self.colors):
                    j=0
            elif isinstance(colors,list):
                if j >= len(colors):
                    j=0
            elif isinstance(colors,str):
                colors = [colors]
                j=0
            if colors == None:
                plt.errorbar(means[i][0],means[i][1],yerr = nstd * np.array(stds[i][1]),color = self.colors[j],ls = linestyle,lw = linewidth,marker = marker,label = labels[i]) # plot with y-axis error bars
            else:
                if clr.is_color_like(colors[j]):
                    plt.errorbar(means[i][0],means[i][1],yerr = nstd*np.array(stds[i][1]),color = colors[j],ls = linestyle,lw = linewidth,marker = marker,label = labels[i])
                else:
                    print("*** WARNING ***: '{0}' is not recognized as a valid color code".format(colors[j]) )
                    plt.errorbar(means[i][0],means[i][1],yerr = nstd * np.array(stds[i][1]),color = self.colors[j],ls = linestyle,lw = linewidth,marker = marker,label = labels[i])
                    colors = None
            j+=1
        if is_legend:
            plt.legend(numpoints=1,frameon=True,loc=legend_location)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)


    def AverageDistributionsCI(self,means,stds,nstd,datatype,labels,colors,title,xlabel,ylabel,is_legend,legend_location):
        """
        Plots the average and standard deviation.

        Input:
         - *means* (nested list)
         - *stds* (nested list)
         - *nstd* (float)
         - *labels* (list)
         - *linestyle* (string)
         - *linewidth* (float)
         - *marker* (string)
         - *colors* (list)
         - *title* (string)
         - *xlabel* (string)
         - *ylabel* (string)
         - *is_legend* (boolean)
         - *legend_location* [default = 'upper right'] (string/integer)
        """
        assert nstd > 0, "Error: The number of STDs must be a value larger than zero"
        plt.figure(self.plotnum)
        datatype_indices = [labels.index(Id) for Id in datatype]
        for i in datatype_indices:
            L_s_amount = copy.copy(means[i][0])
            L_mu =  copy.copy(means[i][1])
            L_sigma =  copy.copy(stds[i][1])

            # Add an additional value
            L_s_amount.append(L_s_amount[-1]+1)
            L_mu.append(L_mu[-1])
            L_sigma.append(L_sigma[-1])

            X_i = []
            Y_i = []
            L_errors = []
            for j in range(len(L_s_amount)):
                if (not L_s_amount[j] == L_s_amount[0]) and (not L_s_amount[j] == L_s_amount[-1]):
                    X_i.append(L_s_amount[j])
                    Y_i.append(L_mu[j-1])
                    L_errors.append(L_sigma[j-1])
                X_i.append(L_s_amount[j])
                Y_i.append(L_mu[j])
                L_errors.append(L_sigma[j])
            X_e = np.concatenate([X_i, X_i[::-1]])
            Y_e = np.concatenate([np.array(Y_i) - nstd*np.array(L_errors) ,(np.array(Y_i) + nstd*np.array(L_errors))[::-1]])

            if colors == None:
                if j >= len(self.colors):
                    j=0
            elif isinstance(colors,list):
                if j >= len(colors):
                    j=0
            elif isinstance(colors,str):
                colors = [colors]
                j=0
            if colors == None:
                plt.fill(X_e-0.5,Y_e,  alpha=.25, ec='None', label='{0} STD confidence interval'.format(nstd),color = self.colors[j])
                plt.plot(np.array(X_i)-0.5,np.array(Y_i),color = self.colors[j])
            else:
                if clr.is_color_like(colors[j]):
                    plt.fill(X_e-0.5,Y_e,  alpha=.25, ec='None', label='{0} STD confidence interval'.format(nstd),color = colors[j])
                    plt.plot(np.array(X_i)-0.5,np.array(Y_i),color = colors[j])
                else:
                    print("*** WARNING ***: '{0}' is not recognized as a valid color code".format(colors[j]) )
                    plt.fill(X_e-0.5,Y_e,  alpha=.25, ec='None', label='{0} STD confidence interval'.format(nstd),color = self.colors[j])
                    plt.plot(np.array(X_i)-0.5,np.array(Y_i),color = self.colors[j])
                    colors = None
        if is_legend:
            plt.legend(numpoints=1,frameon=True,loc=legend_location)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
