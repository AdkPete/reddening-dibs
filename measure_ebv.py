
'''
This script is the one to use to get E(B-V) from a spectrum (or spectra).
Will automatically attempt to fit all of the Friedmann (2011)
DIB features as well as the Na D lines. 

The design is to step through the lines one by one, and allow the user
to decide if the line should be used, and to then set the parameters
for measuring the line EW.

At the end of this, all of the data will be saved to a log file, and the 
derived EW measurements are combined into a E(B-V) measurement with an
uncertainty

'''

import argparse

## Arrays of central wavelenghths.
DIBs = [ ]
Na = []

def measure_spectrum(filename):
    '''
    Function to measure EW for every DIB/Na feature in a spectrum
    '''
    
    return 0

def measure_all_spectra(file_list):
    '''
    Function to measure EW from each spectrum provided
    '''
    
    return 0

def measure_dib_ebv(EW_data):
    '''
    Function to measure E(B-V) from DIB measurements
    '''
    
    return 0

def measure_na_ebv(EW_data):
    
    '''
    Function to measure E(B-V) from Na measurements
    '''
    
    return 0

def measure_ebv(EW_data):
    
    '''
    Function to measure E(B-V) from DIB and Na measurements
    '''
    
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Measure Line EW')
    parser.add_argument('DIB', metavar='DIB_measurement_flag', type=bool, nargs = 1,
                    help='Flag that determines if the script will measure DIB EW')
    
    parser.add_argument('Na', metavar='Na_measurement_flag', type=bool, nargs = 1,
                    help='Flag that determines if the script will measure DIB EW')
    
    parser.add_argument('filenames', metavar='List_of_filenames'
                        , type = str , nargs = '*'
                        , help = "List of filenames that you want to analyze")
    
    args = parser.parse_args()