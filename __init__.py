from . import utils, fitting
import sys

# set options
def SetOptions(calibrations:dict, errors:dict, options:dict, scheme:str):
    
    # set scheme
    options['scheme'] = scheme

    # ranges
    if scheme.lower() == 'sanders':
        options['range'] = [7.0, 8.4]
    elif scheme.lower() == 'bian':
        options['range'] = [7.8, 8.4]
    
    # misc options
    options['verbose'] = False

    # set calibrations and errors
    if scheme.lower() == 'sanders':
        # set calibrations
        calibrations['O3']    = utils.SandersO3
        calibrations['O2']    = utils.SandersO2
        calibrations['R23']   = utils.SandersR23
        calibrations['O32']   = utils.SandersO32
        calibrations['Ne3O2'] = utils.SandersNe3O2
        # set errors
        errors['O3']    = 0.02
        errors['O2']    = 0.07
        errors['R23']   = 0.02
        errors['O32']   = 0.09
        errors['Ne3O2'] = 0.07
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
        errors['O3']    = 0.00
        errors['O2']    = 0.00
        errors['R23']   = 0.00
        errors['O32']   = 0.00
        errors['Ne3O2'] = 0.00
        # conversions
        options['x_to_oh'] = utils.Bian_logOH
        options['oh_to_x'] = utils.Sanders_x
    else:
        print('-> [zcal]: incorrect calibration scheme provided. Exiting...')
        sys.exit()

def Initialise(scheme:str) -> None:
    SetOptions(calibrations, calib_errors, options, scheme)

# set parameters
calibrations = {}
calib_errors = {}
options = {}
scheme = ''
