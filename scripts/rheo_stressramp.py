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

        data_frame = Rstressramp.splitaction_ramp(data_input)

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

        # Find the location of the desired action, and save to a data frame
        # If the action name is not found, it prints an error message
        try: 
            select_data = data_frame.loc[(data_frame[action_header]==iaction[0])]            
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
        # less than 10000%
        ind_nonrupture = np.where(strain < 1e5)[0]
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


        
