"""
Module to process and analyse rheology data 
containing time sweeps / polymerisation curves

Created: March 17th, 2020
Author: Cristina MT
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

class Rtimedep():
    """
    Class with the functions relevant to time sweeps or 
    polymerisation curves.

    Main focus in extracting data from .csv files,
    averaging multiple curves, and fittin the 
    elastic modulus [t]. 
    """

    def readcsv_full2time(filename, 
                        export = True,
                        file_export = None,
                        variables = ['sample des', 'time (sample)', 'modulus', 'force'],
                        sep = ',', dec = '.'):

        """
        Function to select the desired data from raw .csv files.

        INPUT:
            filename :  string containing the file to read
            export :    if True, the selected data is exported to a .csv file
            file_export : name of the file where the data will be exported.
                            if None, then attaches the suffix '_clean_timedep'
                            to the filename.
            variables : array of strings, containing the desired variables.
                        The name of the variables matches - at least partially - 
                        the columns in the imported .csv files. Matching is
                        case insensitive. If more than one column matches a
                        given variable name, all the corresponding columns are
                        included.
            sep :       string indicating the character used as a delimiter in the
                        .csv file. Default is comma ','
            dec :       string indicating the character used for decimal separator
                        in the .csv file. Default is dot '.'

        OUTPUT:
            select_data :   data frame with the selected data. Only returns the value
                        if export  = False. When export = True, the function
                        only exports the data without returning any values.

        """

        # Import the file as a data frame
        data_input = pd.read_csv(filename,
                                sep = sep, decimal = dec)

        print('\n Successfully imported the file:')
        print(filename)

        # Because there is more than one action in the file, 
        # select only the data for time sweep
        # TO DO: make this selection optional for the user

        data_frame = Rtimedep.splitaction_time(data_input)

        # Find the columns that match the desired variable names
        # and select the data within.
        columns = []
        for ivar in variables:
            print('\n Variable to search:', ivar)
            column_names = [x for x in data_frame.columns if ivar in x.lower()]
            print('Variables found:', column_names)
            columns.extend(column_names)

        select_data = data_frame[columns]

        # Export the data to the file specified in file_export or
        # return the data frame if export == False.
        if export == True:
            if file_export == None: 
                file_export = filename.replace('.csv','_clean_timedep.csv')
            select_data.to_csv(file_export, 
                                index=False, 
                                sep = sep, decimal = dec)
            print('\n Selected data exported to:', file_export)
        else:
            return select_data

    def splitaction_time(data_frame, 
                        action_header = 'Action Name', 
                        action_name = 'oscillation'):

        """
        Function to extract the time sweep data from
        a file with multiple types of measurement

        INPUT
            data_frame :    pandas data frame with the full data
            action_header : string with the name of the column
                            containing the type of measurement, or action,
            action_name :   string with the name of the measurement,
                            or action. It accepts a partial match,
                            and is case insensitive.

        OUTPUT
            select_data :   pandas data frame containing only the
                            time sweep data.
        """

        print('\n Splitting data by action name: ', action_name)
        
        # Gets all the actions within the data frame
        iaction = [x for x in data_frame[action_header].unique() if action_name in x.lower()]

        # Find the location of the desired action, and save to a data frame
        # If the action name is not found, it prints an error message
        try: 
            select_data = data_frame.loc[(data_frame[action_header]==iaction[0])]            
        except IndexError: 
            print('ERROR: Action name not found') 
            select_data = None   
        return select_data

    def readtime_sample(data_frame, sample_name,
                        sample_header = 'Sample Description',
                        variables = ['time', 'elastic'],
                        verbose = True):
        """
        Function to extract the temporal variables in the data frame,
        for a given sample. The default is to extract only the 
        time and the elastic modulus.

        INPUT
            data_frame : pandas data frame containing the data
            sample_name : string with the name of the sample
            sample_header : string indicating the name of the column
                            containing the sample name
            variables : array of string with the desired variables.
                        The name can be a partial match to the column 
                        names, and is case insensitive.
            verbose :   If True, it prints the sample name and the
                        column names that match the variables.

        OUTPUT
            extracted_data : numpy array with all the desired variables 
        """

        # Find the columns matching the variable names
        columns = []

        for ivar in variables:
            column_names = [x for x in data_frame.columns if ivar in x.lower()]
            columns.extend(column_names)
        
        if verbose == True:
            print('\n Sample name: ', sample_name)
            print('Extracting variables: ', columns)

        # Extract only the required columns from the data, for the indicated sample
        data_sample = data_frame.loc[(data_frame[sample_header] == sample_name)]
        data_sample = data_sample.dropna()  # Some measurements have a NA value,
                                            # corresponding to the change of 
                                            # measurement type, or action.

        # Store the data in a numpy array and return it.
        extracted_data = []
        for ivar in range(len(variables)):
            extracted_data.append(np.array(data_sample[columns[ivar]].transpose()))

        return extracted_data

    def plot_all(data_frame, t = 'time', yvariable = 'elastic',
                sample_header = 'Sample Description',
                axis_labels = ['Time (s)', 'G\' (Pa)'], show_legend = True,
                save_fig = False, fig_size = (9,5),
                color_scheme = 'viridis',
                linestyle = '-', linewidth = 1.5, marker = None):

        """
        Function to plot the time sweep data for all the samples
        within a file.

        INPUT
            data_frame : pandas data frame containing the full data
            t :         string indicating the name of the time variable
                        to look in the column names.
                        It can be a partial match, case insesitive.
            yvariable : string indicating the name of the dependent
                        variable to look in the column names. It can be
                        a partial match, case insensitive. The default
                        is set to look for the elastic modulus.
            sample_header : String with the column name where the
                        name of the sample is.
            axis_labels : [xlabel, ylabel] containing the names of the
                        labels for the axis, as string.
            show_legend : if True, shows the legend on the plot
            save_fig :  if True, saves the figure. 
                        TO DO: not implemented yet
            fig_size :  (width, height) indicating the size of the figure
            color_scheme : String indicating the color map to be used for
                        plotting the curves. It should match a matplotlib name.
            linestyle : String wiht the linestyle to be used for plotting
                        the curves. It should follow matplotlib.pyplot specs.
            linewidth : float indicating the line width for plotting the curves
            marker :    string indicating the type of marker to use in the plot. 
                        It should match matplotlib.pyplot specs. Default is NONE
            """

        # Read all the sample names in the data frame    
        all_samples = data_frame[sample_header].unique()
        N = len(all_samples);
        
        # Initialize figure, and color scheme
        plt.figure(figsize = fig_size)
        pcm = plt.get_cmap(color_scheme)
        color_lines =  pcm(np.linspace(0,1,N+1))
        ic = 0

        # Loop to plot all the curves, there is a pause of 1 second in between
        # to visualize in 'real time'

        for isample in all_samples:
            print(isample)
            ic = ic+1
            # Extract the time and yvar from the full dataframe, for the 
            # selected sample
            [time, yvar] = Rtimedep.readtime_sample(data_frame, sample_name = isample,
                                            sample_header = sample_header,
                                            variables = [t, yvariable])
            plt.plot(time, yvar, label = isample, c = color_lines[ic],
                    ls = linestyle, lw = linewidth, marker = marker)
            plt.pause(1)
        
        # Add axis labels and legend (if show_legend == True)
        plt.xlabel(axis_labels[0])
        plt.ylabel(axis_labels[1])

        if show_legend == True: plt.legend()

    def mean_polycurve(data_frame, yvariable, tmin, tmax,
                        avg_window = 10, show_plot = True,
                        sample_header = 'Sample Description'):

        """
        Function to average multiple time sweeps. Note that 
        it doesn't use any interpolation, but instead, it averages
        the data points within a time window.

        INPUT
            data_frame : pandas data_frame containing the full data
            yvariable : string with the variable name to be averaged.
                        It can be a partial match, case insensitive.
            tmin : float, the low bound time limit to consider
            tmax : float, the high bound time limit to consider
            avg_window : float, time window to be used for averaging
            show_plot : If True, plots the results
            sample_header : string, column name where the sample name is.

        OUTPUT
            time_mean : numpy array with the average curve of time.
            yvar_mean : numpy array with the average curve of the yvariable
            yvar_std :  numpy array with the standard deviation curve of the yvariable
        """

        # Read all samples within the data frame
        all_samples = data_frame[sample_header].unique()
        N = len(all_samples)

        time_mean = []
        yvar_mean = []
        yvar_std = []

        # Loop to average all the curves within the time window
        for dt in range(tmin, tmax, avg_window):
            time_all = []
            yvar_all = []
            for isample in all_samples:
                # It extracts the time, and yvariable from the data frame
                # for a given sample
                [t, yvar] = Rtimedep.readtime_sample(data_frame, sample_name = isample,
                                                sample_header = sample_header,
                                                variables = ['time', yvariable],
                                                verbose = False)
                # Selects the data within the time window and stores it
                ind_selec = (t>(dt-avg_window)) & (t<=(dt))
                time_all.extend(t[ind_selec])
                yvar_all.extend(yvar[ind_selec])
            # Convert from list to numpy array for mean and isnan to work
            time_all = np.array(time_all)
            yvar_all = np.array(yvar_all)
            try:
                # Get the mean curve, only for non nan values
                # The loop above returns nan values if there's no data
                # in the time window.
                time_mean.append(np.mean(time_all[~np.isnan(time_all)]))
                yvar_mean.append(np.mean(yvar_all[~np.isnan(yvar_all)]))
                yvar_std.append(np.std(yvar_all[~np.isnan(yvar_all)]))
            except TypeError: print('Error in mean calculation')

        # Convert from list to numpy array
        time_mean = np.array(time_mean)
        yvar_mean = np.array(yvar_mean)
        yvar_std = np.array(yvar_std)

        # Plot the average curve and standard deviation
        if show_plot == True:
            Rtimedep.plot_meanpoly(time_mean, yvar_mean, yvar_std)

        return [time_mean, yvar_mean, yvar_std]
        
    def plot_meanpoly(time_mean, yvar_mean, yvar_std, color_scheme = 'gray',
                        axis_labels = ['Time (s)', 'G\' (Pa)'],
                        linestyle = '-', linewidth = 1, alpha_fill = 0.8,
                        fig_size = (9, 5), save_fig = False, marker = None):

        """
        Function to plot the average curve of the time sweep.

        INPUT
            time_mean : numpy array, time (t)
            yvar_mean : numpy array, average curve for yvar[t]
            yvar_std :  numpy array, standard devaition for yvar[t]
            color_scheme : string, color name used as a base for plotting
            axis_labels : list [xlabel, ylabel] with strings containing the
                        axis labels
            linestyle : string, line style to be used for the mean curve.
                        It should match matplotlib.pyplot specs.
            linewidth : float, line width to be used for the mean curve.
                        It should match matplotlib.pyplot specs
            alpha_fill : float, between [0, 1] alpha value (opacity) 
                        to be used for the standard deviation area
            fig_size : tuple, (width, height) for figure size
            save_fig : if True, saves the figure. 
                        TO DO: not implemented yet
            marker : string, marker to be used for the mean curve, 
                        it should match matplotlib.pyplot specs.
        """


        color_fill = 'light' + color_scheme
        color_line = 'dark' + color_scheme
        plt.figure(figsize = fig_size)       
        plt.fill_between(time_mean, yvar_mean - yvar_std, yvar_mean + yvar_std,
                            color = color_fill, alpha = alpha_fill)
        plt.plot(time_mean, yvar_mean, c = color_line,
                            ls = linestyle, lw = linewidth, marker = marker)
        
        plt.xlabel(axis_labels[0])
        plt.ylabel(axis_labels[1])

    def fit_polycurve(time, modulus, exp_fit = 'single', 
                    verbose = True, show_fit = True):

        """
        Function to fit the polymerisation curve with a single
        or double exponential curve, to extract the <Go>

        INPUT
            time : numpy array, time data (t)
            modulus : numpy array, modulus data G'[t] or G''[t]
            exp_fit : 'single' or 'double', indicates the type of fit to be used.
            verbose : if True, prints the results of the fit
            show_fit : if True, plots the results of the fit.

        OUTPUT
            gp_fit : float, Go extracted from fit
            gp_fit_error : float, error in Go
            time_fit : float, tao extracted from fit
            time_error : float, error in tao
        """

        # Initial parameters for fit
        y0 = 0
        A1 = np.max(modulus)
        t1 = 5
        t2 = 2000

        # Convert time to minutes
        time_min = time/60

        # Choose function to fit, single or double exp.
        if exp_fit == 'single': 
            p0 = [y0, A1, t1]
            func_fit = FitFunction.exp_single
        elif exp_fit == 'double':
            p0 = [y0, A1, t1, t2]
            func_fit = FitFunction.exp_double

        # Non linear fit
        pfit, pvar = curve_fit(func_fit, time_min, modulus, p0=p0)
        gfit = func_fit(time_min, *pfit)

        # Extract fit values of interest
        gp_fit = np.round(pfit[0], 3)
        gp_fit_error = np.round(np.sqrt(pvar[0][0]), 3)
        time_fit = np.round(pfit[2], 3)
        time_error = np.round(np.sqrt(pvar[2][2]), 3)

        # Print the results
        if verbose == True: 
            print('\n Fit results summary')
            print('Storage modulus (Pa): ', gp_fit)
            print('Error storage modulus: ', gp_fit_error)
            print('Time (min): ', time_fit)
            print('Error time: ', time_error)

        # Plot the results
        if show_fit == True:
            plt.figure()
            plt.plot(time_min, modulus, c = 'k')
            plt.plot(time_min, gfit, c = 'm')
            plt.xlabel('Time (min)')
            plt.ylabel('G\' (Pa)')
            plt.pause(0.2)

        return [gp_fit, gp_fit_error, time_fit, time_error]

    def fitall_polycurve(file_list, save_results = True, show_fit = True,
                        sample_header = 'Sample Description'):
        """
        Function to fit all the polymerisation curves in multiple files

        INPUT
            file_list :     List of .csv file names containing the data
            save_results :  if True, saves the result to a .txt file
            show_fit :      if True, plots the fitted curve
            sample_header : string, column name where the sample names are.
        OUTPUT
            saves the results of the fit in a .txt file if the option is selected.

        """

        # Name of the file where the results will be saved. 
        # Note that it will be saved in the current working directory
        # TO DO: Generalize the file name

        file_export = 'Results_fit.txt'

        # Initializes the headers of the .txt file
        if save_results == True:
            with open(file_export, 'w') as f:
                f.write('File , Sample , G\' (Pa) , err G\', t0 (min) , t0e')

        # Reads the .csv files within the file_list and fits the data
        for ifile in file_list:
            data_frame = pd.read_csv(ifile)
            # For each sample within the selected file
            for isample in data_frame[sample_header].unique():
               
                # Extracts only the desired data, default now for G'
                # TO DO: Generalise for other yvar
                [t, gp] = Rtimedep.readtime_sample(data_frame, isample, 
                                        verbose = False)

                # Fits the data, and plots the output if selected
                [g0, g0e, t0, t0e] = Rtimedep.fit_polycurve(t, gp, 
                                    show_fit = show_fit)

                # Appends the results in the .txt file, if selected
                if save_results == True:
                    str2write = '\n ' + ifile + ' , ' + isample + ' , ' + \
                                str(g0) + ' , ' + str(g0e) + ' , ' + \
                                str(t0) + ' , ' + str(t0e) + ' , '
                    with open(file_export, 'a') as f:
                        f.write(str2write)

class FitFunction():
    """
    Class with the functions required to fit the time sweep
    or polymeristion curves data from rheology.
    """

    def exp_single(t, *p):
        """
        Function of a single exponential curve of the form
        """
        y0, A1, t1 = p
        y = y0 + A1*np.exp(-t/t1)
        return y  

    def exp_double(t, *p):
        """
        Function of a single exponential curve of the form
        """
        y0, A1, t1, t2 = p
        y = y0 + A1*np.exp(-t/t1) + A1*np.exp(-t/t2)
        
        return y  
