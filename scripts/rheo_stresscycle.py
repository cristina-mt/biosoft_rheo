""""
Module to process and analyse rheology data
containing stress cycles

Note that some basic operations are contained
within the Rstressramp module; here we deal with
operations relative ONLY to stress ramp CYCLES.

Created: April 1st, 2020
Author: Cristina MT
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Rstresscycle():
    """
    Class with the functions relevant to stres cycles
    
    Main focus on computation and data visualization
    (Extracting data is done on Rstressramp)

    TO DO: easier integration with other modules
    """
    
    def plot_kcycles(data_frame, x_header, k_header,
                    cycle_header, win_rolling = 1, **kwargs):

        """
        Function to plot, in log scale, the differential
        elastic modulus, k' for cyclic stress (un)loads
        as a function of stress, strain, or both.

        INPUT
            data_frame :    pandas data frame containing the data
            x_header :      string, name of the column with the x data
            k_header :      string, name of the column with the k data
            cycle_header :  string, name of the column with the cycle labels
            win_rolling :   integer, width of the window for rolling mean
            **kwargs :      plot options for the load cycles and unload cycles

        **kwargs OPTIONS:
            fig_size :     tuple, size of the figure
            N_colors:      integer, number of colors for color scheme
            color_load :   string, color sceheme to be used in load cycles
            lw_load :      float, line width
            ls_load :      string, line style
            marker_load :  string, marker type
            ms_load :      integer, marker size
            mec_load :     string, color scheme for the marker edge. It can
                           take 4 values: 'as_line', 'ligher', 'darker', 'none'
            mfc_load :     string, color scheme for the marker face color.
                           It can take the same values as mec_load

            All of the arguments can be used for unload cycles, by
            substituting '_load' with '_unload'
        """

        # Initialize values for plot options
        try: fig_size = kwargs['fig_size']
        except KeyError: fig_size = (6,6)
        
        load_options = {
                        'color_load': 'viridis', 
                        'lw_load' : 1,
                        'ls_load' : '-',
                        'marker_load' : 'o',
                        'ms_load' : 5,
                        'mec_load' : 'as_line',
                        'mfc_load': 'as_line'
                        }
            
        unload_options = {
                        'color_unload': 'viridis', 
                        'lw_unload' : 1,
                        'ls_unload' : '-',
                        'marker_unload' : 'o',
                        'ms_unload' : 3,
                        'mec_unload' : 'as_line',
                        'mfc_unload': 'none'
                        }

        # Replace default values by user input ones
        for key in kwargs:
            if key in load_options.keys(): 
                load_options[key] = kwargs[key]
            elif key in unload_options.keys():
                unload_options[key] = kwargs[key]
            
        # Initialize figure and set up color scheme
        cycles = data_frame[cycle_header].unique()
        load_cycles = [x for x in cycles if 'Load' in x]
        unload_cycles = [x for x in cycles if 'Unload' in x]
        
        try: N = kwargs['N_colors']
        except KeyError: N = len(load_cycles)

        plt.figure(figsize = fig_size)
        
        pcm_load = plt.get_cmap(load_options['color_load'])
        pcm_unload = plt.get_cmap(unload_options['color_unload'])

        # color scheme for load cycles
        color_lines =  pcm_load(np.linspace(0,1,N+3))
        color_load = color_lines[1:]

        if load_options['mec_load'] == 'as_line':
            color_mec = color_lines[1:]
        elif load_options['mec_load'] == 'lighter':
            color_mec = color_lines[:]
        elif load_options['mec_load']  == 'darker':
            color_mec = color_lines[2:]
        elif load_options['mec_load'] == 'none':
            color_mec = ['None'] * len(color_lines)
        
        if load_options['mfc_load'] == 'as_line':
            color_mfc = color_lines[1:]
        elif load_options['mfc_load'] == 'lighter':
            color_mfc = color_lines[:]
        elif load_options['mfc_load']  == 'darker':
            color_mfc = color_lines[2:]
        elif load_options['mfc_load'] == 'none':
            color_mfc = 'None'
        

        # color scheme for unload cycles
        ucolor_lines =  pcm_unload(np.linspace(0,1,N+2))
        color_unload = ucolor_lines[1:]

        if unload_options['mec_unload'] == 'as_line':
            ucolor_mec = ucolor_lines[1:]
        elif unload_options['mec_unload'] == 'lighter':
            ucolor_mec = ucolor_lines[:]
        elif unload_options['mec_unload']  == 'darker':
            ucolor_mec = ucolor_lines[2:]
        elif unload_options['mec_unload'] == 'none':
            ucolor_mec = 'None'
        
        if unload_options['mfc_unload'] == 'as_line':
            ucolor_mfc = ucolor_lines[1:]
        elif unload_options['mfc_unload'] == 'lighter':
            ucolor_mfc = ucolor_lines[:]
        elif unload_options['mfc_unload']  == 'darker':
            ucolor_mfc = ucolor_lines[2:]
        elif unload_options['mfc_unload'] == 'none':
            ucolor_mfc = ['None'] * len(ucolor_lines)
        
        
        # Iterate over the load cycles
        ic = 0
        for icycle in load_cycles:
            print(icycle)
            ic = ic + 1
            data_cycle = data_frame.loc[data_frame[cycle_header] == icycle]
            data_cycle = data_cycle.rolling(win_rolling).mean()
            xvar = data_cycle[x_header]
            k = data_cycle[k_header]
            plt.plot(xvar, k, c = color_load[ic],
                    lw = load_options['lw_load'],
                    ls = load_options['ls_load'],
                    marker = load_options['marker_load'],
                    ms = load_options['ms_load'],
                    mec = color_mec[ic], 
                    mfc = color_mfc[ic])
            plt.pause(0.5)
        
        # Iterate over the unload cycles
        ic = 0
        for icycle in unload_cycles:
            ic = ic + 1
            print(icycle)
            data_cycle = data_frame.loc[data_frame[cycle_header] == icycle]
            data_cycle = data_cycle.rolling(win_rolling).mean()
            xvar = data_cycle[x_header]
            k = data_cycle[k_header]
            plt.plot(xvar, k, c = color_unload[ic],
                    lw = unload_options['lw_unload'],
                    ls = unload_options['ls_unload'],
                    marker = unload_options['marker_unload'],
                    ms = unload_options['ms_unload'],
                    mec = ucolor_mec[ic], 
                    mfc = ucolor_mfc[ic])
            plt.pause(0.5)

        plt.loglog()
        
        
