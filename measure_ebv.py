
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

import numpy as np

import argparse
import utility as util
import measure_line as ml

## Arrays of central wavelenghths.
DIBs = ['5487.7' , '5705.1' , "5780.5" , "5797.1" , "6196.0" , "6204.5" , "6283.8" , "6613.6"]
Na =[ '5889.98' , '5895.93'] ##D2 , D1
debug = [5780.5]

def measure_spectrum(filename):
    '''
    Function to measure EW for every DIB/Na feature in a spectrum
    '''
    
    spec = ml.spectrum(filename)
    lines = []
    if args.DIB:
        lines += DIBs
    if args.NA:
        lines += Na
    lines = debug
    ewdata = []
    for lc in lines:
        
        lambda_c = float(lc)
        ew , ewsigma = ml.measure_ew(spec , lambda_c , param_file_name = "params.txt")
        print (ew , ewsigma)
        ewdata = [ lc , ew , ewsigma]
        
    return ewdata

def measure_all_spectra(file_list):
    '''
    Function to measure EW from each spectrum provided
    '''
    results = {}
    for fname in file_list:
        ewdata = measure_spectrum(fname)
        results[fname] = ewdata
        
    return results


def build_scatter_table(filename = "data_files/Friedmann_Table_1.txt"):
    
    '''
    Function designed to read in the Friedman (2011) DIB data, then build
    a table with the relevant scatters based on that
    
    Parameters
    ----------
    filename : str
        name of the file containing the Friedman table 1 data
    
    Returns
    -------
    scatter : dict
        Dictionary containing bins with the scatter in the Friedman data
        Dictionary keys are a string, set by the dib wavelength
        Bin info is stored as dib_bin
        
    '''
    

    dib_data = util.read_file(filename , delim="\t" , comment = "#" , min_length = 5)
    dibs = ['5487.7' , '5705.1' , "5780.5" , "5797.1" , "6196.0" , "6204.5" , "6283.8" , "6613.6"]
    dib_fits = []
    
    scatter = {}
    
    dib_res = []
    ebv_res = []
    err_res = []
    
    for dibkey in dibs:
        ebv = []
        dib = []
        dib_err = []
        keys = []
        for i in range(len(dib_data["HD"])):
            
            if dib_data[dibkey][i] == "":
                continue
            if "<" in dib_data[dibkey][i] or ">" in dib_data[dibkey][i]:
                continue
            ebv.append(float(dib_data["E_B-V^a"][i]))
            keys.append(dibkey)
            dib.append(float(dib_data[dibkey][i].split(" ")[0]))

            if len(dib_data[dibkey][i].split(" ")) == 1:
                dib_err.append(0)
            dib_err.append(float(dib_data[dibkey][i].split(" ")[-1]))
        keys = np.array(keys)
        dib = np.array(dib)
        dib_err = np.array(dib_err)
        
        ebv = np.array(ebv)
        
        res = []
        fit_ebv = []
        for i in range(len(keys)):
            res = []
            res.append({})
            
            res[-1]["DIB"] = [keys[i]]
            res[-1]["EW"] = [dib[i] / 1000.0]
            res[-1]["Err"] = [dib[i] / 1000.0] ##Placeholder. Doesn't get used
            fit_ebv.append(measure_ebv(res)[0])
        
        
        fit_ebv = np.array(fit_ebv)
        
        ii = np.argsort(dib)
        dib = dib[ii]
        dib_err = dib_err[ii]
        ebv = ebv[ii]
        fit_ebv = fit_ebv[ii]
        
        ebv -= fit_ebv
        fit_ebv -= fit_ebv
        
        N_per_bin = 30
        i = 0
        
        bins = []
        sigmas = []
        data = []
        while i < len(dib):
            
            bin = []
            if i == 0:
                bin.append(0)
            else:
                bin.append(dib[i])
                
            if i + 2 * N_per_bin >= len(dib):
                end = len(dib)
                bin.append(np.inf)
            else:
                end = i+N_per_bin
                bin.append(dib[end])
            
            bins.append(bin)
            sigmas.append(np.std(ebv[i:end]))
            data.append(dib[i:end])


            if end == len(dib):
                break
            i += N_per_bin
        
        scatter[dibkey] = [ bins , sigmas , data ]
        
    return scatter

def get_scatter(ew , dibkey , scatter_dict):
    
    '''
    Function to find appropriate scatter for a given DIB and EW
    '''
    
    bins , sigmas , data = scatter_dict[dibkey]
    for i in range(len(bins)):
        bin = bins[i]
        if ew >= bin[0] and ew < bin[1]:
            return sigmas[i]
        
def measure_dib_ebv(EW_data):
    '''
    Function to measure E(B-V) from DIB measurements
    '''


    ebv = []
    ebvsigmas = []
    for spectrum_fname in EW_data.keys():
        for line in EW_data[spectrum_fname]:
            
            if line[1] < 0:
                continue
            
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
    
    parser.add_argument('NA', metavar='Na_measurement_flag', type=bool, nargs = 1,
                    help='Flag that determines if the script will measure DIB EW')
    
    parser.add_argument('filenames', metavar='List_of_filenames'
                        , type = str , nargs = '*'
                        , help = "List of filenames that you want to analyze")
    
    args = parser.parse_args()
    
    EW_Data  = measure_all_spectra(args.filenames)
    measure_dib_ebv(EW_Data)