
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
import scipy.signal as signal


class spectrum:

    '''
    Class designed to store a spectrum and perform basic
    calculations on it. Includes the EW width codes,
    and some basic continuum fitting using
    Legendre polynomials
    '''

    def __init__(self ,fname):

        '''
        Initizlize a spectrum given the filename containing the 
        data.
        '''

        wave , spectrum = read_spectrum(fname)
        self.wave = wave
        self.spec = spectrum
        self.fname = fname
        self.cont = np.ones(spectrum.shape)
        self.cnorm = False
        self.ii = np.where(np.isfinite(self.spec))
        self.med_smooth = None
        self.delta_lam = self.wave[1] - self.wave[0]
        self.params = None
        self.cont_ii = np.where(wave > 0)
        self.not_cont = np.where(wave < 0 )
        self.int_cont = None
        self.int_ii = None
        self.delta_lambda = 1.4
        
    def get_range(self , wlow , whigh):
        '''
        Filter the spectrum down to only wavelegths between the given
        values. All internal methods will automatically use the current
        cuts.
        Parameters
        ----------
        wlow : float
            Lower wavelength limit (Angstroms)
        whigh : float
            Upper wavelength limit (Angstroms)

        '''
        self.ii = np.where( (self.wave > wlow) & (self.wave < whigh) )

    def read_parameters(self , fname = "params.txt"):
        '''
        reads in the parameter file for interactive mode
        '''
        
        res = {}
        f = open(fname)
        for i in f.readlines():
            if i[0] == "#":
                continue
            elif i.strip() == "":
                continue

            sl = i.split("\t")
            ## Remove empty strings
            sl = [x for x in sl if x]
            if len(sl) == 1:
                continue
            if "other_exclutions" not in sl[0]:
                res[sl[0].strip()] = sl[1]
            else:
                ex_arr = []
                for i in sl[1].split(","):
                    ex_arr.append( [ float(i.split("-")[0]) , float(i.split("-")[1])] )
                res[sl[0].strip()] = ex_arr
                
            continue
        self.params = res
        
        self.delta_lambda = float(self.params["delta_lambda"])

    def plot_spectrum(self , filename = None , pause = False , vline = None , figsize = None , xlims = None , show_cont = True):

        '''
        Function to plot and display this spectrum
        '''
        
        if pause:
            plt.clf()
        if figsize != None:
            plt.figure(figsize = figsize)
        if xlims is not None:
            plt.xlim(xlims)
        plt.ylim(min(self.spec[self.ii]) * 0.9 , max(self.spec[self.ii]) * 1.1) 
        plt.plot(self.wave[self.ii] , self.spec[self.ii] , color = "blue")
        plt.scatter(self.wave[self.ii][self.not_cont] , self.spec[self.ii][self.not_cont] , color = "red" , s = 5)
        if self.int_cont is not None:
            plt.plot(self.wave[self.int_ii] , self.int_cont[self.int_ii] , color = 'red')
        if (max(self.cont) != 1 or min(self.cont) != 1) and not self.cnorm and show_cont:
            plt.plot(self.wave[self.ii] , self.cont[self.ii] , color = "orange")
        if vline is not None:
            for val in vline:
                plt.axvline(val , color = "red" , ls = "--")
        if self.cnorm:
            plt.plot(self.wave[self.ii] , np.ones(self.cont[self.ii].shape) , color = "orange")
        plt.xlabel("Wavelength ($\AA$)")
        plt.ylabel("Normalized Spectrum")

        if filename is None and not pause:
            res = plt.ginput(2)
            plt.show(block=False)
            vals = [res[0][0] , res[1][0]]
            low = min(vals)
            high = max(vals)
            return low , high
        elif pause:
            plt.pause(0.5)
        else:
            plt.savefig(filename)
            plt.close()
        return 0 , 0
    
    def fit_legendre_continuum(self , N = 8):

        '''
        Function to fit the continumm using an Nth order Legendre polynomial.
        Default is to use 8th order following Friedman (2011)
        saves our fit continumm under self.continuum

        Parameters
        ----------
        N : int
            order of the polynomial to fit

        Returns
        -------
            norm_spec : np.array
                Continuum-normalized spectrum
            
        '''
        median_filtered_spec = signal.medfilt(self.spec , 101)
        
        f1 = np.polynomial.Legendre( [1] + (N-1) * [0] , domain = (min(self.wave[self.ii][self.cont_ii]) , max(self.wave[self.ii][self.cont_ii]) ) )
        res = f1.fit(self.wave[self.ii][self.cont_ii] , median_filtered_spec[self.ii][self.cont_ii] , N)
        
        cont = res(self.wave[self.ii])
        self.cont[self.ii] = cont
        
        self.norm_spec = self.spec[self.ii] / self.cont[self.ii]

        return self.spec[self.ii] / self.cont[self.ii]
    

    def get_continuum_ii(self , center_w , low , high):
        
        ex_width = 100
        if self.params == None:
            print ("Warning : No parameter file provided, defaulting to standard values.")
            print ("These defaults may not make sense, please look at the spectrum")
            exclusions =  [ [ low , high] ]
        else:
            ex_width = float(self.params["width"])
            exclusions = [ [ low , high] ] + self.params['other_exclutions_fit']
        
        wave = self.wave[self.ii]
        condition = (abs(center_w - wave) < ex_width )

        for i in exclusions:
            ii = np.where(condition)
            
            condition[ii] = np.logical_or( ( wave[ii] < i[0] ) , ( wave[ii] > i[1] ) )
        #ii = np.where(condition)
        #condition[ii] = ( self.wave[ii] > high )
        
        self.cont_ii = np.where(condition)
        self.not_cont = np.where( np.logical_not(condition) )
        
        
            
    def measure_ew(self,  center_w , wlow , whigh , use_cont = False , debug = False):

        '''
        Uses direct integration to measure a line's equivalent width

        Errors are from Vollmann and Eversbergh (2006)
        The delta_Lambda methodology is based on https://arxiv.org/pdf/1108.1083 (See appendix A)
        
        Parameters
        ----------
        center_w : float
            central wavelength of the emission line
        wlow : float
            Low point used for continuum integration
        whigh : float
            High point used for continuum integration
        use_cont : bool
            If true, integrate under the continuum. Intended for use with
            a continuum normalized spectrum

        Returns
        -------
        EW : float
            Equivalent width
        sigma_ew : float
            Uncertainty on Equivalent Width
        '''

        
        if use_cont and not self.cnorm:
            raise ValueError("Error : Cannot integrate under cont. without first continuum normalizing the spectrum")


        ii = np.where( ( self.wave >= wlow - (self.delta_lambda / 2.0))
                       & (self.wave <= whigh + (self.delta_lambda / 2.0) ) )
        Fbar = np.mean(self.spec[ii])
        
        if use_cont:
            cont = self.cont[ii]
        else:
            i1 = np.where( ( self.wave <= wlow) & (self.wave >= wlow - self.delta_lambda ) )
            c_low = np.mean(self.spec[i1])
            
            i2 = np.where( ( self.wave >= whigh) & (self.wave <= whigh + self.delta_lambda ) )
            c_high = np.mean(self.spec[i2])
            cont_slope = (c_low - c_high) / (wlow - (self.delta_lambda / 2.0) -  whigh + - (self.delta_lambda / 2.0))
            
            cont_b = c_low - cont_slope * (wlow - (self.delta_lambda / 2.0))
            cont = self.wave * cont_slope + cont_b
            self.int_cont = cont
            self.int_ii = ii

        Fcbar = np.mean(cont[ii])
        integrand = (1 - self.spec[ii] / cont[ii])
        
        
        ### Next up, we have numerical integration of this:

        W = 0
        height = self.delta_lam 
        for i in range(len(integrand)):
            if i == len(integrand) - 1:
                continue
            W += ( integrand[i] + integrand[i+1] ) * height / 2.0
            
        if debug:
            
            print (W)
                
            plt.plot(self.wave[ii] , self.spec[ii])
            plt.plot(self.wave[ii] , cont[ii] , color = 'red')
            plt.xlabel("Wavelength ($\AA$)")
            plt.ylabel("Normalized Flux")
            plt.show()
        
        
        snr = np.mean(self.spec[self.ii][self.cont_ii]) / np.std(self.spec[self.ii][self.cont_ii] - self.cont[self.ii][self.cont_ii])
        
        sigma = np.sqrt(1 + Fcbar / Fbar) *  ( (whigh - wlow) - W ) / snr
        
        return W , sigma
    
    def measure_ew_grid(self , center_w , Nsamp = 100):

        '''
        Function to run on continuum normalized spectrum
        Estimates the EW using a grid of end points for the integration
        Uses the distribution or resulting values to estimate the uncertainty
        on the equivalend with

        Parameters
        ---------
        center_w : float
            centeral wavelength of the line in question
        Nsamp : int
            Number of wavelength combos to try

        '''

        dw = abs( self.wave  - center_w )

        center_i = np.where( dw == np.min(dw))[0][0]
        
        max_lw = center_i - 1
        min_hw = center_i + 1

        if self.spec[center_i] > self.cont[center_i]:
            return -1 , -1 , -1 , -1
        while True:
            if self.spec[max_lw] < self.cont[max_lw]:
                max_lw -= 1
            if self.spec[min_hw] < self.cont[min_hw]:
                min_hw += 1

            if self.spec[min_hw] > self.cont[min_hw] and self.spec[max_lw] > self.cont[max_lw]:
                break
        
        if self.wave[min_hw] - self.wave[max_lw] < 1:
            return -1 , -1 , -1 , -1
        #max_lw -= 1
        #min_hw += 1
        self.plot_spectrum(vline = [center_w , self.wave[min_hw] , self.wave[max_lw]] , pause = False)
        N = int(np.sqrt(Nsamp)) + 1

        base_ew , base_sigma = self.measure_ew( center_w , self.wave[max_lw] , self.wave[min_hw] , debug = False )

        EWs = []
        for i in range(N):
            for k in range(N):
                
                EW , sigma = self.measure_ew( center_w , self.wave[max_lw - i ] , self.wave[min_hw + k ] )
                EWs.append(EW)

        print (np.mean(EWs) * 1000 , np.std(EWs) * 1000 , base_ew * 1000 , base_sigma * 1000)
        #plt.hist(EWs)
        #plt.show()
        return np.mean(EWs) , np.std(EWs) , base_ew , base_sigma
    
    def get_limits(self , center_w):
        if self.params["low"].strip() != "0":
           
            try:
                low = float(self.params["low"])
                high = float(self.params["high"])
                return low , high
            except:
                pass
            
        dw = abs( self.wave  - center_w )
        center_i = np.where( dw == np.min(dw))[0][0]
        
        max_lw = center_i - 1
        min_hw = center_i + 1
        
        if self.spec[center_i] > self.cont[center_i]:
            low = center_w - 2
            high = center_w + 2
        else:
            while True:
                if self.spec[max_lw] < self.cont[max_lw]:
                    max_lw -= 1
                if self.spec[min_hw] < self.cont[min_hw]:
                    min_hw += 1
                if min_hw > len(self.spec) - 2 or max_lw < 0:
                    break
                if self.spec[min_hw] > self.cont[min_hw] and self.spec[max_lw] > self.cont[max_lw]:
                    break
        try:
            low = self.wave[max_lw]
            high = self.wave[min_hw]
        except:
            breakpoint()

        return low , high
    
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
    return wl , flux

def measure_ew(s1 , lambda_c , param_file_name = "params.txt"):
    '''
    Main function that measures an EW with an uncertainty
    Parameters
    __________
    s1 : spectrum object
         object containing the spectrum that you'd like to analyze
    lambda_c : float
        central wavelength for the line (in Angstroms)
    '''


    
    
    if not plt.fignum_exists(1):
        fig = plt.figure(figsize = (8,8))
    else:
        plt.clf()

    ## Read parameter file
    s1.read_parameters( param_file_name )
    
    ## Cut down to within 20 AA of central wavelength
    s1.get_range(lambda_c - 20 , lambda_c + 20)
    
    ## Fit continuum model
    
    s1.get_continuum_ii(lambda_c , lambda_c-3 , lambda_c + 3 )
    
    s1.fit_legendre_continuum(int(s1.params["L_order"]))

    ## Get integration limits
    low , high = s1.get_limits(lambda_c)

    
    mlow = None
    mhigh = None
    while True:
        if mlow is None:
            s1.read_parameters( "params.txt")
        
        s1.get_continuum_ii(lambda_c , low , high)
        s1.fit_legendre_continuum(int(s1.params["L_order"]))
        
        if mlow is None:
            low , high = s1.get_limits(lambda_c)
        EW , Sigma = s1.measure_ew(lambda_c , low , high)
        print ("Equivalent Width = {} +/- {}\n".format(EW , Sigma))
        pw = float(s1.params["Plot_Limits"])
        s1.plot_spectrum(pause = True , vline = [lambda_c , low , high] ,
             xlims = [ lambda_c - pw , lambda_c + pw ] , show_cont = True)
        
        u = input("Use this measurement? (y/n/r) --> ")
        if u == "": ##Repeat calculation, including param file read
            continue
        elif u.lower()[0] == "n":
            low , high = s1.plot_spectrum(vline = [lambda_c , low , high] , 
            xlims = [ lambda_c - pw , lambda_c + pw ] , show_cont = True)
            mlow = low
            mhigh = high
            
        elif u.lower()[0] == "y":
            return EW , Sigma
        
        elif u.lower().strip()[0] == "r":
            mlow = None
            mhigh = None
  
        
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Measure Line EW')
    parser.add_argument('fname', metavar='Filename', type=str, nargs = 1,
                    help='filename for your spectrum')
    parser.add_argument('lambda_c', metavar='lambda', type=float, nargs=1,
                    help='central wavelength for your line')
    args = parser.parse_args()

    spec1 = spectrum(args.fname[0])
    
    EW , sigma_EW = measure_ew(spec1, args.lambda_c[0])
    print (EW , "+/-" , sigma_EW)
    