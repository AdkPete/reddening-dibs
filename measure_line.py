
'''
Main function for measuring equivalent widths and their uncertainties.
Does not include any of the conversions from EW to E(B-V)
You can use this by running this script with the name of the spectrum file
as the first argument, followed by the central wavelength of the line that you
would like to measure.
'''


import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits

def read_ascii_spectrum(fname):
    '''
    Simple function to read SMARTS spectra out of the provided ASCII files.
    Assumed formatting is an ASCII file containing two columns, lambda first,
    then the flux. Columns should be seperated by spaces (any #).
    
    Parameters
    ----------
    fname : str
        string containing the name of the ascii file to read
    
    Returns
    --------
    wls : arr
        array containing the wavelengths of each point in the spectrum
    flux : arr
        array containg the fluxes for each point in the spectrum
    '''
    
    f = open(fname)
    
    wls = []
    flux = []
    for i in f.readlines():
        if i[0] == "#":
            continue
        sl = i.split(" ")
        wla = False
        for k in sl:
            if k == "":
                continue
            if not wla:
                wls.append(float(k))
                wla = True
            else:
                flux.append(float(k))
                break
            
    wls = np.array(wls)
    flux = np.array(flux)
    
    return wls , flux
    
def read_fits_spectrum(file):
    
    '''
    simple function to read in a spectrum
    for now, official support is limited to the ARAS specra,
    although many other spectra utilize similar FITS files
    
    Parameters
    ----------
    filename : str
        path to the file containing the desired spectrum

    Returns
    -------
        (wlarr , spectrum)
        wlarr : np.array
            Array containing wavelengths associated with spectrum in Angtroms
        spectrum : np.array
            Array containing our actual spectrum, normalized by the average
            flux. Has no units due to this normaliztion.

    '''
    

    hdu = fits.open(file)
    spectrum = hdu[0].data
    
    wl0 = hdu[0].header['CRVAL1']
    delta_wl = hdu[0].header['CDELT1']

    
    wlarr = np.array(range(len(spectrum))).astype(np.float64)
    wlarr *= delta_wl
    wlarr += wl0

    
    return wlarr , spectrum / np.mean(spectrum)

def read_spectrum(filename):
    '''
    Function to read in a spectrum and auto-detect file type
    Will call the correct function and then return a wavelength array
    and a flux array, which may be used for equivalent width measurements.
    
    Parameters
    __________
    filename : string
        name of spectrum file
    Returns
    _______
        wlarr : array of floats
            wavelengths for the spectrum (in Angstroms)
        flux : array of floats
            flux values (arbitrary units, normalized to mean flux)
    '''
    if ".fit" in filename:
        wl , flux = read_fits_spectrum(filename)
    else:
        try:
            wl , flux = read_ascii_spectrum(filename)
        except:
            raise ValueError("Failed to read in spectrum")
def measure_ew(wlength , flux , lambda_c):
    '''
    Main function that measures an EW with an uncertainty
    Parameters
    __________
    wlength : array of floats
         wavelength array for the spectrum (in Angstroms)
    flux : array of floats
        flux array for the spectrum (arbitrary units)
    lambda_c : float
        central wavelength for the line (in Angstroms)
    '''
    return 0

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Measure Line EW')
    parser.add_argument('fname', metavar='Filename', type=str, nargs = 1,
                    help='filename for your spectrum')
    parser.add_argument('lambda_c', metavar='lambda', type=float, nargs='+',
                    help='central wavelength for your line')
    args = parser.parse_args()

