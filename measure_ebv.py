
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
import os

import argparse
import utility as util
import measure_line as ml
from tabulate import tabulate
## Arrays of central wavelenghths.
DIBs = ['5487.7', '5705.1', "5780.5", "5797.1", "6196.0", "6204.5",
        "6283.8", "6613.6"]
Na = [ '5889.98' , '5895.93'] ##D2 , D1
debug = [5780.5]


### Model parameters from Friedmann 2011 / Poznanski 2012

fit_par = {}
fit_par["5487.7"] = [-6.41E-02 , 9.672E-03]
fit_par["5705.1"] = [-1.74E-01,1.2E-02]
fit_par["5780.5"] = [-8.36E-03 , 1.98E-03]
fit_par["5797.1"] = [-2.86E-02 , 5.74E-03]
fit_par["6196.0"] = [-5.07E-02 , 2.11E-02]
fit_par["6204.5"] = [-7.22E-02 , 5.99E-03]
fit_par["6283.8"] = [ -7.71E-02 , 9.57E-04]
fit_par["6613.6"] = [1.96E-02 , 4.63E-03]
fit_par["5895.93"] = [-1.76 , 2.47] ##5895.93
fit_par["5889.98"] = [-1.91 , 2.16] ##5889.98

unc_par = {} ##Parameter uncertainties
unc_par["5487.7"] = [1.31E-02 , 0.25E-03]
unc_par["5705.1"] = [0.16E-01, 0.16E-02]
unc_par["5780.5"] = [3.48E-03 , 0.01E-03]
unc_par["5797.1"] = [0.57E-02 , 0.06E-03]
unc_par["6196.0"] = [0.56E-02 , 0.06E-02]
unc_par["6204.5"] = [0.67E-02 , 0.08E-03]
unc_par["6283.8"] = [ 0.78E-02 , 0.17E-04]
unc_par["6613.6"] = [0.37E-02 , 0.04E-03]
unc_par["5895.93"] = [0.17 , 0]
unc_par["5889.98"] = [0.15 , 0]

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
    #lines = debug
    ewdata = []
    for lc in lines:
        
        lambda_c = float(lc)
        ew , ewsigma = ml.measure_ew(spec , lambda_c , param_file_name = "params.txt")
        print (ew , ewsigma)
        
        ewdata.append([ lc , ew , ewsigma])
        
        
    return ewdata

def measure_all_spectra(file_list):
    '''
    Function to measure EW from each spectrum provided
    '''
    if len(file_list) == 1 and os.path.isdir(file_list[0]):
        files = os.listdir(file_list[0])
        for i in range(len(files)):
            files[i] = os.path.join(file_list[0] , files[i])
        
    else:
        files = file_list
    results = {}
    for fname in files:
        try:
            ewdata = measure_spectrum(fname)
        except:
            continue
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
    

    dib_data = util.read_file(filename, delim="\t", comment = "#",
        min_length = 5)
    dibs = ['5487.7', '5705.1', "5780.5", "5797.1", "6196.0", "6204.5",
        "6283.8", "6613.6"]

    scatter = {}

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
            res = {}

            res["f1"] = [[]]
            res["f1"][0].append(keys[i])
            res["f1"][0].append(dib[i] / 1000.0)
            res["f1"][0].append(dib[i] / 1000.0)
            
            fit_ebv.append(measure_ebv(res,use_scatter=False)[0])
        
        
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
        

def EW_to_ebv(EW , Err , lc):
    '''
    Function to convert EW to E(B-V) using the Friedmann 2011
    calibrations (Poznanski 2012 in the case of the Na lines)
    '''

    fit_key = str(lc)

    if fit_key not in fit_par.keys():
        print (f"Error, line {fit_key} is not supported\
            (not in fit parameters)")
        raise ValueError

    if fit_key in DIBs:

        a = fit_par[fit_key][0]
        b = fit_par[fit_key][1]
        sa = unc_par[fit_key][0]
        sb = unc_par[fit_key][1]
        bv = a + b * EW * 1000.0
        measure_err = Err * 1000.0 * b
        return bv , measure_err

    elif fit_key in Na:
        a = fit_par[fit_key][0]
        b = fit_par[fit_key][1]
        sa = unc_par[fit_key][0]
        sb = unc_par[fit_key][1]
        bv = 10 ** (a + b * EW)
        
        measure_err = b * np.log(10) * bv * Err * 1000.0
        
        return bv , 0.35 * bv
    
    else:
        print (f"Error, line {fit_key} is not supported \
            (not a known DIB or Na Line)")
        raise ValueError

def measure_ebv(EW_data , use_scatter = True , useDIB = True , useNa = True):
    '''
    Function to measure E(B-V) from DIB measurements
    '''
    
    if use_scatter and useDIB:
        ## If not DIB, then no need to build the DIB scatter table
        scatter = build_scatter_table()

    ebv = []
    ebv_err = []
    table_rows = []
    header = ["Line" , "EW (milliangstroms)" , "EW Error (milliangstroms)",
        "E(B-V)" , "E(B-V) Error"]
    table_rows.append(header)
    lines = []
    if useDIB:
        lines += DIBs
    if useNa:
        lines += Na
        
    if len(lines) == 0:
        print ("Error, no lines to measure")
        return 0

    for lc in lines:
        EW_values = []
        EW_errors = []
        
        for specname in EW_data.keys():
            for line in EW_data[specname]:
                
                
                if line[1] < 0: ## Negative EW measurements ; not useful here
                    continue
                if str(line[0]) == lc:
                    EW_values.append(line[1])
                    EW_errors.append(line[2])
                    
        if len(EW_values) == 0: ## No EW Measurements, continue
            continue 
        
        
        
        EW_values = np.array(EW_values)
        EW_errors = np.array(EW_errors)
        ## Convert to a E(B-V) Value. Starts with mean EW

        w = 1 / np.array(EW_errors) ** 2
        EW = ( np.sum(w * EW_values) / np.sum(w) )
        
        ## W/ Standard Error prop
        Err = 1 / np.sqrt( np.sum( w ) )
        
        bv , measure_err = EW_to_ebv(EW , Err , lc)
        ebv.append(bv)
        
        
        if not use_scatter:
            ebv_err.append(1)
            continue
        
        if str(lc) in DIBs:
            sigma_ebv = get_scatter(EW * 1000.0 , str(lc) , scatter)
        
            joint_err = np.sqrt(sigma_ebv ** 2 + measure_err ** 2 )
        elif str(lc) in Na:
            joint_err = np.sqrt(measure_err ** 2 + (bv * 0.35) ** 2)
        else:
            print (f"Warning: Line {lc} not known")
        
        ebv_err.append(joint_err)
        
        row = [] ## Row for nice table output
        if lc not in Na:
            row.append(lc)
        else:
            if lc == Na[0]:
                row.append(lc + " Na D2")
            else:
                row.append(lc + " Na D1")
        row.append(round(EW * 1000.0 , 1))
        row.append(round( measure_err * 1000.0 , 1))
        row.append(round(bv , 2))
        row.append(round(joint_err , 2))
        table_rows.append(row)
    ### Now compute the final E(B-V)
    ebv_arr = np.array(ebv)
    ebv_errors = np.array(ebv_err)
    
    w = 1/ ebv_errors ** 2
    ebv = np.sum(w * ebv_arr) / np.sum(w)
    
    ## W/ Standard Error prop + some algebra,
    ebv_err = 1 / np.sqrt( np.sum( w ) )
    
    return ebv , ebv_err , ebv_arr , ebv_errors , table_rows



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Measure Line EW')
    parser.add_argument('DIB', metavar='DIB_measurement_flag', type=bool,
        nargs = 1,
        help='Flag that determines if the script will measure DIB EW')
    
    parser.add_argument('NA', metavar='Na_measurement_flag', type=bool,
        nargs = 1,
        help='Flag that determines if the script will measure DIB EW')
    
    parser.add_argument('filenames', metavar='List_of_filenames_or_directory',
        type = str , nargs = '*', help = "List of filenames (or a directory\
        containing these files) that you want to analyze")
    
    args = parser.parse_args()
    
    EW_Data  = measure_all_spectra(args.filenames)
    
    EBV, Err, ebv1, ebv2, table = measure_ebv(EW_Data, useDIB = args.DIB,
        useNa = args.NA)
    
    print (f"\n\nE(B-V) = {round(EBV , 2)} +/- {round(Err , 2)}")
    
    print (tabulate(table , tablefmt="fancy_grid"))