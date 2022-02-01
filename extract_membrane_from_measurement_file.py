# coding=utf-8
# Read and analyse pulse measurements.
# (c) Oisin Creaner, DIAS, 24th January, 2022

# based on extract membrane hits from measurement.pro by Gerhard Ulbricht, January 2014
# based on  pulsefit.pro by Ben Mazin, June, 2013

# imports the argparse library for parsing command line arguments
import argparse

# imports the json library for reading and writing json files
import json

# imports the warnings library
import warnings

# imports the OS library for disc operations
import os


def arg_parser():
    """
    This function parses the arguments from the command line and returns a dictionary containing information used by the
    rest of the program to do its job

    Several options are provided: Positional arguments, followed by optional
    arguments followed by interactive entry of the argument values.

    future expansions to arguments will allow the user to specify modes of
    operation and the type of output generated
    """

    # creates an argparse variable that the user can input arguments to at the command line
    parser = argparse.ArgumentParser()

    # specifies the arguments

    ###############################################################################
    # Input masterfile
    ###############################################################################

    # adds an optional argument for the masterfile
    parser.add_argument("--masterfile", "-M", default=None,
                        help='''
    Path to the Masterfile as used by the original program.  This will likely be dropped for future versions  
                        ''')

    # collects all of the arguments from parser now that they have been given their values
    args = parser.parse_args()

    # creates an output directory to store the required information to control the program that will be used throughout
    # the rest of the system
    modes = {}

    # extracts the parsed arguments from args to modes
    modes['masterfile']=args.masterfile

    return (modes)

def read_json(filename):
    """
    This function reads in a json file and returns its content

    :param filename: a path to a file
    :return: jsonfile_content, a dictionary containing the information contained in jsonfile
    """
    # check if the file exists
    if os.path.exists(filename):
        # if it does, open it
        with open(filename) as jsonfile:
            print('reading in jsonfile '+filename)
            # and read it
            jsonfile_content = json.load(jsonfile)
    # Otherwise the file does not exist
    else:
        warnings.warn("Specified file "+filename+" does not exist, continuing with blank dictionary")
        # return an empty dictionary which can fit into later function calls on it
        jsonfile_content = {}


    return (jsonfile_content)

def fit_pulses(modes):
    """
    procedure to call the least square fit for every single pulse
    :return:
    """
    
    warnings.warn('fit_pulses currently empty.')
    pass


def generate_iqsweepfit():
    """
    fit sweep data = resonator IQ loop, but only if the fit file doesn't already exist
    :return:
    """
    ####################################################################################################################
    # From the comments in the annotated version, it looks like this isn't the way we're going to do this going forward#
    # I'll hold off on developing it further until after further discussion.  If you find this comment in a live version
    # this conversation probably ended in a decision that this would never be needed.                                  #
    ####################################################################################################################
    warnings.warn('generate_iqsweepfit currently empty.')
    pass
def main ():
    # main body of the program.  As much as possible, this will call other functions

    ######################################
    # Parsing the command line arguments #
    ######################################
    # runs the argument parser function
    modes = arg_parser()

    ###############################
    # Reading from the masterfile #
    ###############################
    # if the masterfile is not present in the arguments at all
    if 'masterfile' not in modes:
        # return an empty dictionary which can fit into later function calls on it
        warnings.warn('Masterfile not specified in arguments, proceeding without it')
        masterfile_content = {}
    # Otherwise the masterfile argument exists in the input arguments
    else:
        # if the argument is not blank
        if modes['masterfile'] is not None:
            # retrieve the content from the masterfile
            masterfile_content = read_json(modes['masterfile'])

            ####################################################################
            # there probably needs to be a parse masterfile function call here #
            ####################################################################

        # otherwise the argument exists and is not blank
        else:
            # return an empty dictionary which can fit into later function calls on it
            warnings.warn('Masterfile not specified in arguments, proceeding without it')
            masterfile_content = {}

    ###############################
    # Reading in the iqsweep file #
    ###############################
    # if the iq_sweep_fit has been specified in the arguments
    if 'iq_sweep_fit' in modes:
        # read the contents of the iq_sweep_fit from the file specified in the arguments
        iq_sweep_fit_content = read_json(modes['iq_sweep_fit'])
    # Otherwise, if the iq_sweep_fit has been specified in the Masterfile instead
    elif 'iq_sweep_fit' in masterfile_content:
        # read the contents of the iq_sweep_fit from the file specified in the Masterfile
        iq_sweep_fit_content = read_json(masterfile_content['iq_sweep_fit'])
    # otherwise, return an empty dictionary which can fit into later function calls on it
    else:
        iq_sweep_fit_content = {}

    ######################################
    # Generating IQ sweep data if needed #
    ######################################
    # An empty dictionary evaluates to false, a dictionary with content evaluates to true
    # if the dictionary has content
    if iq_sweep_fit_content:
        # the sweep does not need to be generated
        # warn the user that it was there
        warnings.warn('IQ-loop fit already present')
    # otherwise (the dictionary is empty)
    else:
        # calls functions to generate it
        generate_iqsweepfit(modes)
        # Eoin has a python verion of this to use

    ######################
    # Fitting the pulses #
    ######################
    # begins the process of fitting the pulses
    print ('fitting pulses')
    fit_pulses()

    # report success
    print ('done!')
    pass

# Call the main script if in the base envionment to allow imports
if __name__ == '__main__':

    main()