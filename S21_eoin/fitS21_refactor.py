#note, works for IQ files saved using fsweep macro saved on the VNA
#import usual packages
import warnings

import numpy as np
from scipy.signal import argrelextrema
from scipy.signal import argrelmin
from scipy import optimize
from scipy.optimize import curve_fit
# import pylab
# import matplotlib.pyplot as plt
import array as arr
import csv
#from fpdf import FPDF
# from matplotlib.backends.backend_pdf import PdfPages
from datetime import date
import os
import pandas
from datetime import datetime
import json

def read_S21_data(filename):
    '''
    This file reads in S21 data from a CSV file with headers or JSON file with pulses as appropriate
    :return:
    '''
    # creates an emtpy dictionary that will be returned
    S21_data={}

    extension = os.path.splitext(filename)[1]

    if extension == '.csv':
        S21_data=csv_reader(filename)
    elif extension == '.json':
        S21_data = json_reader(filename)
    else:
        warnings.warn('File type not known: '+extension)

    return (S21_data)

def csv_reader(filename):

    # sets up variables to hold the names of the columns in the IQ data according to the old structure
    freq_name_old='Freq(Hz)'
    I_name_old='S21(REAL)'
    Q_name_old='S21(IMAG)'
    # sets up the number of lines in the header from the old structure
    header_lines=6
    footer_lines=1

    # Extracts the needed information from the filename
    basename,timestamp,timestamp_format=parse_csv_filename(filename)

    # creates a dictionary which will hold the data from the file
    # uses the basename of the file as the name of the data
    data = {'name':basename}

    # creates an empty list of pulses that will allow for graceful failures if there is no suitable data
    data['pulses'] = []

    # Reads the csv file based on the structure defined above
    pulse_data=pandas.read_csv(filename,
                               header=header_lines,
                               skipfooter=footer_lines,
                               engine='python')

    # pulls the data from the appropriate columns
    freq_data = pulse_data[freq_name_old].tolist()
    I_data=pulse_data[I_name_old].tolist()
    Q_data=pulse_data[Q_name_old].tolist()

    # puts that data into the dictionary of pulses
    data['pulses'].append(
        {
        'pulseID': basename,
        'timestamp': timestamp,
        'timestamp_format': timestamp_format,
        'Frequency': freq_data,
        'I': I_data,
        'Q': Q_data
        })

    return data

def json_reader(filename):
    with open(filename) as json_file:
        data = json.load(json_file)
    return data

def parse_csv_filename(filename):
    '''
    This function takes a filename as an argument and returns the base filename and the timestamp extracted from it for
    use in other partsof the program.  As a default, if the filename isn't structured as expected, it returns the
    current date and time instead and raises a warning
    :param filename:
    :return:
    '''

    # this is the easy bit: extract the base filename
    basename=os.path.splitext(filename)[0]

    # Tries to parse out the sweep date from the filename
    try:
        # removes the frequency number from the basename
        sweep_date=basename.rsplit('_',1)[0]

        # parses the date from that
        date_parsed = datetime.strptime(sweep_date, 'sweep_%Y_%m_%d')

        # creates a string for the timestamp and the timestamp format to allow them to be read
        timestamp_format = '%Y%m%d'
        timestamp = date_parsed.strftime(timestamp_format)

    # if that doesn't work: e.g. file name is not formatted correctly
    except ValueError:
        # raises a warning to the user
        warnings.warn('Date not able to be parsed from the filename: '+filename+'\nDefaulting to current time')
        # Creates a timestamp and format string for the current date and time
        timestamp_format = '%Y%m%d%H%M%S'
        timestamp = datetime.now().strftime(timestamp_format)

    #returns what has been parsed
    return basename,timestamp,timestamp_format

def calculate_derived_data(pulse):
    '''
    This fuction takes a pulse (a dictionary containing the I and Q data) and calculates the derived values
    Amplitude, Phase for each measurement (Frequency).

    Values are returned by adding them to the dictionary
    :param pulse:
    :return:
    '''

    # calls the function to calculate the amplitude
    amplitudes=calculate_amplitude(pulse)
    # normalizes magnitude data relative to max value
    normalized_amplitudes = amplitudes / np.max(amplitudes)
    # calculates square magnitude
    square_magnitudes = normalized_amplitudes ** 2
    # converts the amplitudes to a list and adds it to the pulse dictionary
    pulse['Amplitude']=amplitudes.tolist()
    pulse['Normalised Amplitude']=normalized_amplitudes.tolist()
    pulse['Square Amplitude']=square_magnitudes.tolist()

    # calls the function to calculate the phase
    phases=calculate_phases(pulse)
    # converts the phase to a list and adds it to te pulse dictionary
    pulse['Phase']=phases.tolist()

def calculate_amplitude(pulse):
    '''
    This function takes a pulse  (a dictionary containing I and Q data)
    and calculates the amplitudes of each value by using pythagoras theorem
    :param pulse:
    :return:
    '''

    # extracts the I values from the list in the pulse into a numpy array
    I_values = np.array(pulse['I'])

    # extracts the Q values from the list in the pulse into a numpy array
    Q_values = np.array(pulse['Q'])

    # calculates the amplitude from Pythagoras formula applied to I and Q
    amplitudes = np.sqrt((I_values) ** 2 + (Q_values) ** 2)

    # returns the amplitudes
    return amplitudes

def calculate_phases(pulse):
    '''
    This function takes a pulse  (a dictionary containing the I and Q data)
    and calculates the phases using the arctan2 operation in numpy
    :param pulse:
    :return:
    '''

    # extracts the I values from the list in the pulse
    I_values = pulse['I']
    # extracts the Q values from the list in the pulse
    Q_values = pulse['Q']

    # calculates the phases using arctan2 operations via numpy
    phases = np.arctan2(Q_values, I_values)
    return phases

def fit_curve(S21_data):
    return{}

def plot_data(S21_data, fit_data):
    pass


def main():
    filename='sweep_2021_11_10_3871350000.json'
    S21_data=read_S21_data(filename)

    for pulse in S21_data['pulses']:
        calculate_derived_data(pulse)

        fit_data=fit_curve(pulse)

        plot_data(pulse, fit_data)



# Call the main script if in the base envionment to allow imports
if __name__ == '__main__':

    main()