# handle imports
import numpy as np
import matplotlib.pyplot as plt
import zcal
import sys

# -=-=-=- SANDERS 2023 CALIBRATIONS -=-=-=-
# link: 

def SandersO3(x:float) -> float:
    return 0.834 - (0.072 * x) - (0.453 * x * x)

def SandersO2(x:float) -> float:
    return 0.067 + (1.069 * x) 

def SandersR23(x:float) -> float:
    return 1.1017 + (0.026 * x) - (0.331 * x * x)

def SandersO32(x:float) -> float:
    return 0.723 - (1.153 * x) 

def SandersNe3O2(x:float) -> float:
    return -0.386 - (0.998 * x) 

# -=-=-=- BIAN CALIBRATIONS -=-=-=- 
# link: 

def BianO3(x:float) -> float:
    return 43.9836  + (-21.6211 * x) + (3.4277 * x ** 2) + (-0.1747 * x ** 3)
    

def BianO2(x:float) -> float:
    print('-> [ERROR]: Bian O2 calibration does not exist. Returning nan.')
    return np.nan

def BianR23(x:float) -> float:
    return 138.0430 + (54.8284 * x) + (7.2954 * x ** 2) - (0.32293 * x ** 3)

def BianO32(x:float) -> float:
    return (8.54 - x) / 0.59

def BianNe3O2(x:float) -> float:
    return (7.80 - x) / 0.63

# -=-=-=- FITTING PARAMETER CONVERSIONS -=-=-=-
def Bian_logOH(x:np.ndarray) -> np.ndarray:
    return x - 12.

def Bian_x(logOH:np.ndarray) -> np.ndarray:
    return 12. + logOH

def Sanders_logOH(x:np.ndarray) -> np.ndarray:
    return x - 12. + 8.

def Sanders_x(logOH:np.ndarray) -> np.ndarray:
    return 12. + logOH - 8.

# -=-=-=- RUNTIME SETTING METHODS -=-=-=-
def SetCalibration(source:str) -> None:
    zcal.Initialise(source.lower())

# -=-=-=- Dust Attenuation Methods -=-=-=- 

def Cardelli_Attenuation(wl:np.ndarray, EBV:float, Rv:float=3.1) -> np.ndarray:
    """ Calculates the Cardelli Extinction curve, for given wl in μm """
    # convert wl into 1/μm
    x = 1. / wl
    # create Aλ array
    AλAv = np.empty_like(x)

    # infrared
    if min(x) < 0.11: # have to check
        mask = (x >= 0.3) & (x <= 1.1)
        a_x =  0.564 * (np.power(x[mask], 1.61))
        b_x = -0.527 * (np.power(x[mask], 1.61))
        AλAv[mask] = a_x + (b_x / Rv)

    # optical
    mask = (x >= 1.1) & (x <= 3.3)
    y = x[mask] - 1.82 # modified parameter
    a_x = 1. + (0.17699 * y) - (0.50447 * y ** 2) - (0.02427 * y ** 3) + \
          (0.72085 * y ** 4) + (0.01979 * y ** 5) - (0.7753 * y ** 6) + \
          (0.32999 * y ** 7)
    b_x = (1.41338 * y) + (2.28305 * y ** 2) + (1.07233 * y ** 3) - \
          (5.38434 * y ** 4) - (0.62251 * y ** 5) + (5.3026 * y ** 6) - \
          (2.09002 * y ** 7)
    AλAv[mask] = a_x + (b_x / Rv)

    # uv + fuv
    mask = (x >= 3.3) & (x < 5.9)
    a_x = 1.72 - (0.316 * x[mask]) - (0.104 / ((x[mask] - 4.67) ** 2 + 0.341))
    b_x = -3.09 + (1.825 * x[mask]) + (1.206 / ((x[mask] - 4.62) ** 2 + 0.263))
    AλAv[mask] = a_x + (b_x / Rv)
    
    # uv bump
    mask = (x >= 5.9) & (x <= 8.0)
    a_x = 1.72 - (0.316 * x[mask]) - (0.104 / ((x[mask] - 4.67) ** 2 + 0.341)) - \
           (0.04473 * (x[mask] - 5.9) ** 2) - (0.009779 * (x[mask] - 5.9) ** 3)
    b_x = -3.09 + (1.825 * x[mask]) + (1.206 / ((x[mask] - 4.62) ** 2 + 0.263)) + \
           (0.2130 * (x[mask] - 5.9) ** 2) + (0.1207 * (x[mask] - 5.9) ** 3)
    AλAv[mask] = a_x + (b_x / Rv)

    # fuv extension (up to 0.1μm)
    mask = (x >= 8.0)
    y = x[mask] - 8.
    a_x = -1.073 - (0.628 * y) + (0.137 * y ** 2) - (0.070 * y ** 3)
    b_x = 13.670 + (4.257 * y) - (0.420 * y ** 2) + (0.374 * y ** 3)
    AλAv[mask] = a_x + (b_x / Rv)

    # return kλ
    return (AλAv * Rv) * EBV

def DustCorrect(wl:float, wl_grid:np.ndarray, Aλ:np.ndarray) -> float:
    """ Generates the dust correction factor for a given emission line flux """
    return np.power(10, 0.4 * np.interp(x=wl, xp=wl_grid, fp=Aλ))
