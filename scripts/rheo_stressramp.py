"""
Module to process and analyse rheology data
containing stress ramps

Created: March 24th, 2020
Author: Cristina MT
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Rstressramp():
    """
    Class with the functions relevant to stress ramps

    Main focus on extracting data from .csv files,
    computing K' and data visualization
    """

    def readcsv_full2ramp(filename,
                        export = True,
                        file_export = None,
                        action_name = 'viscometry ramp',
                        variables = ['sample des', 'stress', 'strain (sample)'],
                        sep = ',',
                        dec = '.'):

        """
        Function to select the desired data from raw .csv files.
        TO DO: combine this and Rtimedep.readcsv_full2time to a general function

        INPUT
            filename : string, file to read
            export : if True, the selected data is exported to a .csv file
            file_export : string, name of the file where the data will be exported.
                        if None, then attaches the suffix '_clean_stress_ramp'
                        to the file name.
            action_name : string, name of the dataset where the ramp data is
            variables : list of strings, desired variables to be extracted. The 
                        name can be a partial match of the column name, and is 
                        case insensitive. If more than one column matches a given
                        variable name, all the corresponding columns are included.
            sep :       string, character used as a delimiter in the .csv file
            dec :       string, character used as a decimal separator in the .csv file

        OUTPUT
            select_data : data frame with the selected data. Only returns the value
                        if export = False. When export = True, the function only 
                        exports the data without returning any values.
        """

        # Import the file as a data frame
        data_input = pd.read_csv(filename, 
                                sep = sep, decimal = dec)

        print('\n Successfully imported the file: ')
        print(filename)

        # Because there is more than one action in the file,
        # select only the data for the stress ramp
        # TO DO: make this selection optional for the user

        data_frame = Rstressramp.splitaction_ramp(data_input,
                                                  action_name = action_name)

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
                file_export = filename.replace('.csv','_clean_stress_ramp.csv')
            select_data.to_csv(file_export, 
                                index=False, 
                                sep = sep, decimal = dec)
            print('\n Selected data exported to:', file_export)
        else:
            return select_data

    def splitaction_ramp(data_frame, 
                        action_header = 'Action Name', 
                        action_name = 'viscometry ramp'):

        """
        Function to extract the stress ramp data from
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
                            stress ram[] data.
        """

        print('\n Splitting data by action name: ', action_name)
        
        # Gets all the actions within the data frame
        iaction = [x for x in data_frame[action_header].unique() if action_name in x.lower()]

        print(iaction)
        data_frame.set_index(action_header, inplace = True)
        
        # Find the location of the desired action, and save to a data frame
        # If the action name is not found, it prints an error message
        try:
            select_data = data_frame.loc[iaction]
            select_data.reset_index(inplace = True)
        except IndexError: 
            print('ERROR: Action name not found') 
            select_data = None   
        return select_data

    def compute_k(stress, strain, show = None,
                remove_neg = True):

        """
        Function to compute the differential storage modulus
        from the slope of the stress vs strain curve.

        INPUT
            stress : numpy array or list, Shear Stress (in Pa) data
            strain : numpy array or list, Shear Strain (in %) data
            show : 'stress', 'strain', 'both', or None. Plots the result
            remove_neg : if True, removes data where strain is negative
        OUTPUT
            stress : numpy array, mean value of the stress (in Pa) where k is computed
            strain : numpy array, mean value of strain (in %) where k is computed
            k : numpy array, differential storage modulus, (in Pa)
        """

        # Work with numpy arrays
        stress = np.array(stress)
        strain = np.array(strain)

        # Start by cleaning the data from any NaN value
        ind_nan = np.isnan(strain) | np.isnan(stress)
        stress = stress[~ind_nan]
        strain = strain[~ind_nan]

        # Clean the data from values after rupture, strain must be
        # less than 5000%
        ind_nonrupture = np.where(strain < 5e3)[0]
        stress = stress[ind_nonrupture]
        strain = strain[ind_nonrupture]

        # Remove data where strain is negative. Note that if recording
        # the absolute strain of the sample, strain can be negative 
        # in the initial interval. This data is tipically not useful 
        # and therefore not desired.

        if remove_neg == True: 
            ind_positive = np.where(strain >= 0)
            stress = stress[ind_positive]
            strain = strain[ind_positive]

        # Compute the differential values of strain and stress
        diff_stress = stress[1:] - stress[:-1]
        diff_strain = strain[1:] - strain[:-1]

        # Compute k' and the mean values of stress and strain
        k = diff_stress / diff_strain * 100  # multiplied by 100, because strain is in %
        stress = (stress[1:] + stress[:-1])/2
        strain = (strain[1:] + strain[:-1])/2

        # Show the results if desired
        if show == 'stress': Rstressramp.plot_k([stress], k)
        elif show == 'strain': Rstressramp.plot_k([strain], k)
        elif show == 'both': Rstressramp.plot_k([stress, strain], k)
        elif show is not None: print('Error: cannot plot: ', show)

        return [stress, strain, k]

    def plot_k(x, k, linewidth = 1.5, 
                marker = 'o', color = 'k', marker_facecolor = 'k'):

        """
        Function to plot, in log scale, the differential storage modulus, k
        as a function of stress, strain, or both. 

        INPUT
            x : list of numpy arrays of dependent variables
            k : numpy array, differential storage modulus
            linewidth : float, width of the line to plot
            marker : string, marker of the lineplot, needs to be compatible with
                    matplotlib.pyplot
            color : color for the lineplot, and marker border, needs to
                    be compatible with matplotlib.pyplot
            marker_facecolor : color of the marker, compatible with 
                    matplotlib.pyplot
        """

        # Plot the first variable
        x1 = x[0]
        plt.figure(figsize = (9,5))
        plt.plot(x1, k, c = color, lw = linewidth, marker = marker, 
                mec = color, mfc = marker_facecolor)
        plt.loglog()
        plt.ylabel('$K\'$ (Pa)')

        # If there is more than one dependent variable, 
        # Plot also the second variable in a different figure
        try: 
            x2 = x[1]
            plt.xlabel('$\sigma$ (Pa)')
            plt.pause(0.1)
            plt.figure(figsize =(9, 5))
            plt.plot(x2, k, c = color, lw = linewidth, marker = marker, 
                mec = color, mfc = marker_facecolor)
            plt.loglog()
            plt.ylabel('$K\'$ (Pa)')
            plt.xlabel('$\gamma$ (%)')
        except IndexError: pass

    def export_kall(data_frame, file_export = None,
                    remove_neg = True,
                    group_header = 'Sample Description',
                    stress_header = 'Shear stress(Pa)',
                    strain_header = 'Shear strain (sample)(%)'):

        """
        Function to compute the differential storage modulus
        for all the data groups (e.g. samples, interals, experiments)
        within a data_frame

        INPUT
            data_frame : pandas data frame with the full data
            file_export : string, name of the file where data will be exported
                        if None, it saves to 'All_k_curves.csv'
            remove_neg : if True, removes data where strain is negative
            group_header : string, name of the column where the data group label are
            stress_header : string, name of the column where the stress data is
            strain_header : string, name of the column where the strain data is
        OUTPUT
            all_data : data frame with the computed stress, strain, k'
            It also saves the data_rame to file_export.
        """

        groups_all = []
        s_all = []
        y_all = []
        k_all = []

        for igroup in data_frame[group_header].unique():
            data_group = data_frame.loc[data_frame[group_header] == igroup]
            stress = np.array(data_group[stress_header])
            strain = np.array(data_group[strain_header])

            [s, y, k] = Rstressramp.compute_k(stress, strain, remove_neg = remove_neg)

            groups_all.extend([igroup]*len(s))
            s_all.extend(s)
            y_all.extend(y)
            k_all.extend(k)
            
        all_data = pd.DataFrame()
        all_data[group_header] = groups_all
        all_data['Stress (Pa)'] = s_all
        all_data['Strain (%)'] = y_all
        all_data['K prime (Pa)'] = k_all
        
        if file_export is None: file_export = 'All_k_curves.csv'
        all_data.to_csv(file_export, index = False)

        return all_data

    def mean_kall_interp(filename, xvariable,num_interp = 100, show_plot = True,
                        sample_header = 'Sample Description',
                        stress_header = 'Stress (Pa)',
                        strain_header = 'Strain (%)',
                        k_header = 'K prime (Pa)',
                        sep = ',', dec = '.'):

        """
        Function to compute the mean curve for the 
        differential elastic modulus for all the data within a file
        Note that it is based on interpolation!

        INPUT
            filename : string, name of the file with the whole data
            xvariable : string, can be 'stress' or 'strain', indicating 
                    over which variable to compute the mean.
            show_plot : if True, shows the results in a plot
            sample_header : string, name of the column with the sample label is
            stress_header : string, name of the column with the shear data
            strain_header : string, name of the column with the strain data
            sep : string, character used as delimiter in csv file
            dec : string, character used as decimal separator in csv file

        OUTPUT
            xinterp : numpy array, vector used for interpolation
            kmean : numpy array, mean curve of k
            kstd : numpy array, standard deviation curve of k
        """         

         # Read data and get all the samples within the data frame
        data = pd.read_csv(filename, sep = sep, decimal = dec)
        all_samples = data[sample_header].unique()

        # Define which dependent variable to extract
        if 'stress' in xvariable: xvar = stress_header
        elif 'strain' in xvariable: xvar = strain_header

        # Loop to get mean values of minimum and maximum xdata for the samples
        xmin = []; xmax = []
        for isample in all_samples:
            data_sample = data.loc[data[sample_header] == isample]
            xsample = np.array(data_sample[xvar])
            xmin.append(np.min(xsample))
            xmax.append(np.max(xsample))

        xmin_avg = np.mean(np.array(xmin))
        xmax_avg = np.mean(np.array(xmax))
        xmax_std = np.std(np.array(xmax))

        print('Rupture: ', xmax_avg, '+/-', xmax_std)
        # Build interpolation vector
        xmin_log = np.log10(xmin_avg)
        xmax_log = np.log10(xmax_avg)
        xinterp = np.logspace(xmin_log, xmax_log, num = num_interp)

        #Loop to get the interpolated curves for each sample within the file
        k_all = []
        for isample in all_samples:
            data_sample = data.loc[data[sample_header] == isample]
            xsample = data_sample[xvar]
            ksample = data_sample[k_header]
            k_interp = np.interp(xinterp, xsample, ksample)
            k_all.append(k_interp)
            
        k_all = np.array(k_all)
        kmean = np.mean(k_all, axis = 0)
        kstd = np.std(k_all, axis = 0)

        # Plot the average curve and standard deviation, if desired
        if show_plot == True:
            plt.fill_between(xinterp, kmean - kstd, kmean + kstd, color = 'lightgray',
                        alpha = 0.8)
            plt.plot(xinterp, kmean, c = 'darkgray', marker = 'o', mfc = 'w')
            plt.ylabel('$K\'$ (Pa)')
            plt.xlabel(xvar)
            plt.loglog()

        return [xinterp, kmean, kstd]
    


    def mean_kall_window(filename, xvariable, xmin_log = -1, xmax_log = 5, winavg_number = 50,
                    show_plot = True,
                    sample_header = 'Sample Description',
                    stress_header = 'Stress (Pa)',
                    strain_header = 'Strain (%)',
                    k_header = 'K prime (Pa)',
                    sep = ',', dec = '.'):

        """
        Function to compute the mean curve for the 
        differential elastic modulus for all the data within a file
        Note that it is based on window averaging, and not interpolation!

        INPUT
            filename : string, name of the file with the whole data
            xvariable : string, can be 'stress' or 'strain', indicating 
                    over which variable to compute the mean.
            xmin_log : float, minimum value for average -> 10**xmin
            xmax_log : float, minimum value for average -> 10**xmax
            winavg_number : number of windows used to average, in logspace
            show_plot : if True, shows the results in a plot
            sample_header : string, name of the column with the sample label is
            stress_header : string, name of the column with the shear data
            strain_header : string, name of the column with the strain data
            sep : string, character used as delimiter in csv file
            dec : string, character used as decimal separator in csv file

        OUTPUT
            xmean : numpy array, mean value of the xvariable
            kmean : numpy array, mean curve of k
            kstd : numpy array, standard deviation curve of k
        """

        # Read data and get all the samples within the data frame
        data = pd.read_csv(filename, sep = sep, decimal = dec)
        all_samples = data[sample_header].unique()

        # Define which dependent variable to extract
        if 'stress' in xvariable: xvar = stress_header
        elif 'strain' in xvariable: xvar = strain_header

        xmean = []
        kmean = []
        kstd = []

        # Loop to average all the curves within the window

        avg_windows = np.logspace(xmin_log, xmax_log, num = winavg_number)
        avg_windows = [round(x, 3) for x in avg_windows]

        for dw in range(len(avg_windows)-1):
            x_all = []
            k_all = []
            for isample in all_samples:
                # It extracts the xvariable and the k data from the data 
                # frame for a given sample
                data_sample = data.loc[data[sample_header]==isample]
                xdata = data_sample[xvar]
                kdata = data_sample[k_header]
                #Selects the data within the avg window and stores it
                ind_selec = (xdata > avg_windows[dw]) & (xdata <= avg_windows[dw+1])
                x_all.extend(xdata[ind_selec])
                k_all.extend(kdata[ind_selec])
            # Convert list to numpy array for mean and isnan to work properly
            x_all = np.array(x_all)
            k_all = np.array(k_all)
            try:
                # Get the mean curve, only for non values
                xmean.append(np.mean(x_all[~np.isnan(x_all)]))
                kmean.append(np.mean(k_all[~np.isnan(k_all)]))
                kstd.append(np.std(k_all[~np.isnan(k_all)]))
            except TypeError: print('Error in mean calculation')

        # Convert from list to numpy array
        xmean = np.array(xmean)
        kmean = np.array(kmean)
        kstd = np.array(kstd)

        # Plot the average curve and standard deviation, if desired
        if show_plot == True:
            plt.fill_between(xmean, kmean - kstd, kmean + kstd, color = 'lightgray',
                        alpha = 0.8)
            plt.plot(xmean, kmean, c = 'darkgray', marker = 'o', mfc = 'w')
            plt.ylabel('$K\'$ (Pa)')
            plt.xlabel(xvar)
            plt.loglog()

        return [xmean, kmean, kstd]


