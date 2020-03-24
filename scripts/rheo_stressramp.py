"""
Module to process and analyse rheology data
containing stress ramps

Created: March 24th, 2020
Author: Cristina MT
"""

import pandas as pd

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

        
        
