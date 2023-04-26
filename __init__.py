from . import utils, fitting
import sys

# set options
def SetOptions(calibrations:dict, errors:dict, options:dict, wavelengths:dict, scheme:str):
    
    # set scheme
    options['scheme'] = scheme

    # ranges (0.01 to 2 solar)
    if scheme.lower() == 'sanders':
        options['range'] = [6.7, 9.0]
    elif scheme.lower() == 'bian':
        options['range'] = [6.7, 9.0]
    elif scheme.lower() == 'nakajima':
        options['range'] = [6.7, 9.0] # varies for diagnostic

    # misc options
    options['verbose']      = False
    options['corner_plots'] = True
    options['save_pkl']     = True
    options['dust_correct'] = False

    # set all applicable line wavelengths (microns)
    wavelengths['oiii5007'] = 0.500824
    wavelengths['oii'] = (0.372709 + 0.372988)/2.
    wavelengths['hbeta'] = 0.486269
    wavelengths['neiii'] = 0.386981
    
    #wavelengths['oiii4960'] = 0.496030                  

    # set calibrations and errors
    if scheme.lower() == 'sanders':
        # set calibrations
        calibrations['O3']    = utils.SandersO3
        calibrations['O2']    = utils.SandersO2
        calibrations['R23']   = utils.SandersR23
        calibrations['O32']   = utils.SandersO32
        calibrations['Ne3O2'] = utils.SandersNe3O2
        # set errors
        errors['O3']    = 0.09 #0.02
        errors['O2']    = 0.22 #0.07
        errors['R23']   = 0.06 #0.02
        errors['O32']   = 0.29 #0.09
        errors['Ne3O2'] = 0.24 #0.07
        # conversions
        options['x_to_oh'] = utils.Sanders_logOH
        options['oh_to_x'] = utils.Sanders_x
    elif scheme.lower() == 'bian':
        # set calibrations
        calibrations['O3']    = utils.BianO3
        calibrations['O2']    = utils.BianO2 
        calibrations['R23']   = utils.BianR23
        calibrations['O32']   = utils.BianO32
        calibrations['Ne3O2'] = utils.BianNe3O2
        # set errors
        errors['O3']    = 0.10
        errors['O2']    = 0.13
        errors['R23']   = 0.08
        errors['O32']   = 0.19
        errors['Ne3O2'] = 0.20
        # conversions
        options['x_to_oh'] = utils.Bian_logOH
        options['oh_to_x'] = utils.Bian_x
    elif scheme.lower() == 'nakajima':
        # set calibrations
        calibrations['O3']    = utils.NakajimaO3
        calibrations['O2']    = utils.NakajimaO2 
        calibrations['R23']   = utils.NakajimaR23
        calibrations['O32']   = utils.NakajimaO32
        calibrations['Ne3O2'] = utils.NakajimaNe3O2
        # set errors
        errors['O3']    = 0.16
        errors['O2']    = 0.27
        errors['R23']   = 0.10
        errors['O32']   = 0.39
        errors['Ne3O2'] = 0.39
        # conversions
        options['x_to_oh'] = utils.Nakajima_logOH
        options['oh_to_x'] = utils.Nakajima_x
    else:
        print('-> [zcal]: incorrect calibration scheme provided. Exiting...')
        sys.exit()

def Initialise(scheme:str) -> None:
    SetOptions(calibrations, calib_errors, options, wavelengths, scheme)

# set parameters
calibrations = {}
calib_errors = {}
wavelengths = {}
options = {}
scheme = ''
