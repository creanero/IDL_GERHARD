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

import matplotlib.pyplot as plt

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
    # converts the square magnitudes to dB
    square_magnitudesdb = 10 * np.log10(square_magnitudes)
    # converts the amplitudes to a list and adds it to the pulse dictionary
    pulse['Amplitude']=amplitudes.tolist()
    pulse['Normalised Amplitude']=normalized_amplitudes.tolist()
    pulse['Square Magnitude']=square_magnitudes.tolist()
    pulse['Square Magnitude dB']=square_magnitudesdb.tolist()


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

def fit_curve(pulse):
    fit_data={}
    input_parameters=set_parameters(pulse)
    fit_data['popt'], fit_data['popc'] =curve_fit(amplitudeequation,
                                                  pulse['Frequency'],
                                                  pulse['Square Magnitude'],
                                                  input_parameters['p0'],
                                                  sigma = input_parameters['Weighting'],
                                                  absolute_sigma=False,
                                                  maxfev=10000)

    curve_parameters=extract_curve_parameters(fit_data['popt'], fit_data['popc'], pulse)
    fit_data.update(curve_parameters)

    fitted_curve=create_fitted_data(fit_data,pulse)
    fit_data.update(fitted_curve)

    align_curve(pulse)

    return fit_data

def set_parameters(pulse):
    input_parameters={}
    input_parameters['p0']=np.array([1, 1, 1, 1, 10000, np.mean(pulse['Frequency'])]) #initial guesses for fit
    input_parameters['Weighting'] = np.ones(len(pulse['Frequency']))
    return input_parameters

def extract_curve_parameters(popt, pcov, pulse):
    curve_parameters={}
    Q = popt[4]  # total Q is taken from the fit
    Q = abs(Q)
    curve_parameters['Q']=Q

    fr = popt[5]  # resonant frequency is taken from the fit
    # fr_MHz = fr / 10**3 #resonant frequency in MHz
    curve_parameters['fr'] = fr

    Qi = Q / np.min(pulse['Normalised Amplitude'])  # calculates intrinsic Q using equation 28 from Zmuidzinas review
    Qi = abs(Qi)
    curve_parameters['Qi'] = Qi


    Qc = Q / (1 - np.min(pulse['Normalised Amplitude']))  # calculates coupling Q using equation 29 from Zmuidzinas review
    Qc = abs(Qc)
    curve_parameters['Qc'] = Qc

    return curve_parameters

def create_fitted_data(fit_data,pulse):
    popt=fit_data['popt']
    fitted_curve={}

    fitted_curve['Square Magnitude'] = amplitudeequation(pulse['Frequency'],
                                                         popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])  # saves square magnitude fit data

    fitted_curve['Normalised Amplitude'] = np.abs(np.sqrt(fitted_curve['Square Magnitude']))  # takes sqrt to get amplitude

    fitted_curve['Amplitude'] = fitted_curve['Normalised Amplitude'] * np.max(pulse['Amplitude'])  # unnormalizes data relative to max value

    return fitted_curve

def align_curve(pulse):
    centre=get_centre(pulse)
    cable_values=calculate_cable_delay(pulse)
    translated_data=translate_data(pulse, centre, cable_values)
    rotated_data=rotate_data(translated_data)
    fit_phase()

    pass

def get_centre(pulse):
    # get centre of loop
    # find centre of loop
    centre={}
    centre['x'] = (np.max(pulse['I']) + np.min(pulse['I'])) / 2
    centre['y'] = (np.max(pulse['Q']) + np.min(pulse['Q'])) / 2
    return(centre)

def calculate_cable_delay(pulse):
    # Next step, remove cable delay term from data.
    # To do this, normalize by cable delay loop
    cable_values={}
    number_of_values=len(pulse['I'])

    # first find total angle, between first and last point
    # find initial phase
    initial_phase = np.arctan2(pulse['Q'][0], pulse['I'][0])
    # finds the final phase
    final_phase = np.arctan2(pulse['Q'][-1], pulse['I'][-1])

    # creates a linear approximation of the drift between initial and final phase
    # Using linspace prevents floating point errors from accumulating
    # thanks Andrew Jaffe (2009) and BRian Burns (2020) via StackOverflow
    cable_values['Phase']=list(np.linspace(initial_phase,final_phase,number_of_values))


    cable_values['Amplitude'] = np.max(pulse['Amplitude'])  # takes max value
    cable_values['I'] = cable_values['Amplitude'] * np.cos(cable_values['Phase'])  # finds I component of cable
    cable_values['Q'] = cable_values['Amplitude'] * np.sin(cable_values['Phase'])  # finds Q component of cable

    return cable_values

def translate_data(pulse, centre, cable_values):
    translated_data={}

    # subtract cable delay from IQ data to normalize
    translated_data['Irebased'] = pulse['I'] - cable_values['I']
    translated_data['Qrebased'] = pulse['Q'] - cable_values['Q']

    amplitudes_normalized = np.sqrt(translated_data['Irebased'] ** 2 + translated_data['Qrebased'] ** 2)

    # next have to move loop to be centred on (0, 0)

    translated_data['xc'] = (np.max(translated_data['Irebased']) + np.min(translated_data['Irebased'])) / 2
    translated_data['yc'] = (np.max(translated_data['Qrebased']) + np.min(translated_data['Qrebased'])) / 2



    # finds I component of cable
    translated_data['translatedcircleI'] = translated_data['Irebased'] - translated_data['xc']
    # finds Q component of cable
    translated_data['translatedcircleQ'] = translated_data['Qrebased'] - translated_data['xc']

    # Calculates amplitudes by pythagoras
    translated_data['amplitudes_translated'] = np.sqrt(translated_data['translatedcircleI'] ** 2 +
                                                       translated_data['translatedcircleQ'] ** 2)

    return translated_data

def rotate_data(translated_data):
    rotated_data={}

    # find alpha - the argument of zc, the centre of the normalized circle
    alpha = np.arctan2(translated_data['yc'], translated_data['xc'])

    # now rotate everything by alpha
    theta = np.arctan2(translated_data['translatedcircleQ'], translated_data['translatedcircleI'])
    # need to rotate every point in this circle by -alpha
    theta_rotated = theta - alpha  # rotates by alpha
    rotated_data['I_rotated'] = translated_data['amplitudes_translated']  * np.cos(theta_rotated)  # finds I component of cable
    rotated_data['Q_rotated'] = translated_data['amplitudes_translated']  * np.sin(theta_rotated)  # finds Q component of cable

    # amplitudes_rotated = np.zeros(number_of_values)
    rotated_data['amplitudes_rotated'] = np.sqrt(rotated_data['I_rotated']  ** 2 + rotated_data['Q_rotated']  ** 2)

    # plt.scatter(rotated_data['I_rotated'] , rotated_data['Q_rotated'] , s = 5, color = 'y')
    # plt.show()

    rotated_data['thetafinal'] = np.unwrap(theta_rotated)
    return rotated_data

def fit_phase():
    pass

def plot_data(S21_data, fit_data):
    pass

#E.17 from w thesis - will be used later to  amplitude data fit curve
def amplitudeequation(f, A1, A2, A3, A4, Qr, fr):
    return A1 + A2*(f - fr) + ((A3 + A4*(f - fr)) / (1 + 4 * Qr**2  * ((f - fr)/fr)**2 ))

#E.11 from Gao thesis - can be used to fit phase data
def phaseequation(f, theta0, Qr, fr):
    return -theta0 + 2*np.arctan(2*Qr*(1 - (f/fr)))

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