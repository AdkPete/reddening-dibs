import numpy as np
import matplotlib.pyplot as plt
import argparse
from measure_line import spectrum
import scipy.optimize as opt
import scipy.interpolate as interp
import scipy.special as special

def linear_continuum(spec , wlow , whigh):
    '''
    Function used in model fitting.
    Returns the linear continuum model between two endpoints.
    Unlike other cases, this does not extend the edges by delta_lambda
    '''
    
    ii = np.where( ( spec.wave >= wlow - 5 )
                    & (spec.wave <= whigh + 5))

    finterp = interp.interp1d(spec.wave[ii] , spec.spec[ii])
    ii = np.where( ( spec.wave >= wlow )
                    & (spec.wave <= whigh))
    flow = finterp(wlow)
    fhigh = finterp(whigh)
    
    m = (fhigh - flow) / (whigh - wlow)
    b = fhigh - m * whigh
    
    wavels = spec.wave[ii]
    rflux = spec.spec[ii]
    fluxes = wavels * m + b
    def linmod(wl):
        return wl * m + b
    return wavels , fluxes , rflux , linmod
    
def fit_model( spec , model , x0 , low , high , freeze_centers = False , display = False):
    

    if args.Model.lower()[0] == "v":
        npar = 4
    else:
        npar = 3
    
    spec.get_range(low - 10 , high + 10)
    spec.spec /= np.mean(spec.spec[spec.ii])
    
    cont_wave , cont_flux , rflux , linmod = linear_continuum(spec , low , high)
    
    spec.get_range(low - 10 , high + 10)
    x0[2] = -1 / 5 * np.mean(spec.spec[spec.ii])
    model_flux = model(cont_wave , x0)
    
    def f(x):
        if freeze_centers:
            rx = []
            i = 0
            k = 0
            while i < len(x0):
                
                rx.append(x0[i])
                rx.append(x[k])
                rx.append(x[k+1])
                if npar == 4:
                    rx.append(x[k+2])
                i += npar
                k += npar - 1
            
            model_flux = model(cont_wave , rx)
        else:
            model_flux = model(cont_wave , x)
        
        return np.sum( (model_flux + cont_flux - rflux) ** 2 )
    
    if not freeze_centers:
        res = opt.basinhopping(f , x0 , niter = 500, minimizer_kwargs = {"options": {"maxiter" : 5000} })
    else:
        i = 0
        nx0 = []
        while i < len(x0):
            
            nx0.append( x0[i+1] )
            nx0.append( x0[i+2] )
            if npar > 3:
                nx0.append( x0[i+3])
            i += npar
        res = opt.basinhopping(f , nx0 , niter = 500, minimizer_kwargs = {"options": {"maxiter" : 5000} })
    if freeze_centers:
        rx = []
        i = 0
        k = 0
        while i < len(x0):
            rx.append(x0[i])
            rx.append(res.x[k])
            rx.append(res.x[k+1])
            if npar == 4:
                rx.append(res.x[k+2])
            i += npar
            k += npar-1
    else:
        rx = res.x
    model_flux = model(cont_wave , rx)
    if display:
        print (res)
        plt.plot(spec.wave[spec.ii] , spec.spec[spec.ii])
        plt.plot( cont_wave , cont_flux , color = "orange")
        plt.plot( cont_wave , cont_flux + model_flux , ls = "--" , color = "red")
    i = 0
    flux = 0
    
    
    while i < len(rx):
        mx = rx[i:i+npar]
        i += npar
        mflux = model(cont_wave , mx)
        if display:
            plt.plot( cont_wave , cont_flux + mflux , ls = ":" , color = "black")

    if display:
        plt.title("Model Fitting")
        

        plt.show()
    return rx , linmod


def model_ew(s1 , lambda_c , x0 , model , param_file_name = "params.txt"):
    

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

    model_x , linmod = fit_model(s1 , model , 
                x0 , low , high , freeze_centers=True)
    model_x , linmod = fit_model(s1 , model , 
                model_x , low , high , display = True)
    
    if args.Model.lower()[0] == "v":
        npar = 4
    else:
        npar = 3
    EW = 0
    i = 0
    while i < len(model_x):
        def fmod(wl):
            return model(wl , model_x[i:i+npar]) + linmod(wl)
        mew = evaluate_model_ew(low , high , linmod , fmod)
        print (mew)
        EW += mew
        i += npar
    return EW


    
def gauss_model(wl , x):
    '''
    Short function to define a model for gaussian absorption lines
    can manage multiple absorption features, which will be detected
    based on the length of the parameter array
    '''
    
    if len(x) % 3 != 0:
        raise ValueError("Number of parameters in Gauss fit does not make sense")
    i = 0
    flux = 0
    while i < len(x):
        mean = x[i+0]
        sigma = x[i+1]
        norm = x[i+2]
        if norm > 0:
            norm *= -1
        normalization = norm * 1 / (np.sqrt(2 * np.pi * sigma ** 2) )
        flux += normalization * np.exp(-0.5 * ( wl - mean) ** 2 / (sigma ** 2 ) )
        i += 3
        
    return flux

def lorentzian_model(wl , x):
    '''
    Short function to define a model for lorenzian absorption lines
    can manage multiple absorption features, which will be detected
    based on the length of the parameter array
    '''
    
    if len(x) % 3 != 0:
        raise ValueError("Number of parameters in Lorentzian fit does not make sense")
    i = 0
    flux = 0
    while i < len(x):
        x0 = x[i+0]
        gamma = x[i+1]
        norm = x[i+2]
        if norm > 0:
            norm *= -1
        if gamma < 0:
            gamma *= -1
        normalization = norm / (np.pi * gamma )
        flux += normalization / (1 + ( ( wl - x0) / gamma )  ** 2 )
        i += 3
        
    return flux


def voight_model(wl , x):
    '''
    Short function to define a model for lorenzian absorption lines
    can manage multiple absorption features, which will be detected
    based on the length of the parameter array
    '''
    
    if len(x) % 4 != 0:
        raise ValueError("Number of parameters in Voigt fit does not make sense")
    i = 0
    flux = 0
    while i < len(x):
        x0 = x[i+0]
        sigma = x[i+1]
        gamma = x[i+2]
        norm = x[i+3]
        
        
        if norm > 0:
            norm *= -1
        if gamma < 0:
            gamma *= -1
        if sigma < 0:
            sigma *= -1
        res = special.voigt_profile(wl - x0 , sigma , gamma)
        flux += norm * res
        i += 4
        
    return flux

def evaluate_model_ew(low , high , continuum_model , flux_model , N=10000):
    '''
    Function to compute equivalent widths of a flux model compared
    against a linear continuum
    '''
    wls = np.linspace(low , high , N)
    flux_c = continuum_model(wls)
    flux_mod = flux_model(wls)
    
    delta_lam = wls[1] - wls[0]
    integrand = 1 - (flux_mod / flux_c)
    ew = np.sum(integrand * delta_lam)
    return ew

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Measure Line EW')
    parser.add_argument('fname', metavar='Filename', type=str, nargs = 1,
                    help='filename for your spectrum')
    parser.add_argument('lambda_c', metavar='lambda', type=float, nargs=1,
                    help='central wavelength for your line')
    parser.add_argument('N' , metavar='N_models', type=int , nargs = 1 , 
                        help = "Number of model components to fit")
    parser.add_argument('Model' , metavar='Model Type', type=str , nargs='?' ,
                        default =  "Gaussian" , help = "Type of model")
    args = parser.parse_args()

    
    spec1 = spectrum(args.fname[0])
    
    spec1.get_range(args.lambda_c[0] - 3 , args.lambda_c[0] + 3)

    wl_guesses = spec1.plot_spectrum(vline = [args.lambda_c[0]] , 
        xlims = (args.lambda_c[0] - 3 , args.lambda_c[0] + 3) , input_vals=args.N[0])

    plt.close()
    x0 = []
    
    for i in range(args.N[0]):
        x0.append(wl_guesses[i])
        if args.Model.lower()[0] == "v":
            x0.append(0.1)
        
            x0.append(0.1)
        else:
            x0.append(0.3)
        x0.append(-1)
        
    
    if args.Model.lower()[0] == "g":
        model = gauss_model
    elif args.Model.lower()[0] == "l":
        model = lorentzian_model
    elif args.Model.lower()[0] == "v":
        model = voight_model
    else:
        raise ValueError("Model type not regognized, please select from 'g'aussian, 'l'orentzian, or 'v'oight.")
    EW = model_ew(spec1 , args.lambda_c[0] , x0 , model)
    
    print ("Equivalent Width = {}".format(EW))