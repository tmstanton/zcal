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
    return 43.9836  - (21.6211 * x) + (3.4277 * x ** 2) - (0.1747 * x ** 3)
    

def BianO2(x:float) -> float:
    print('-> [ERROR]: Bian O2 calibration does not exist. Returning nan.')
    return -999.0

def BianR23(x:float) -> float:
    return 138.0430 + (54.8284 * x) + (7.2954 * x ** 2) - (0.32293 * x ** 3)

def BianO32(x:float) -> float:
    return (8.54 / 0.59) - (x / 0.59)

def BianNe3O2(x:float) -> float:
    return np.power(10, (-7.80 / 0.63) - (x / 0.63))

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

