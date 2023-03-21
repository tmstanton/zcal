# imports
import zcal
import numpy as np
import scipy.stats as stats
import dynesty
from dynesty import utils as dyfunc

# main method
def FitZg(scheme:str, ratios:np.ndarray, errors:np.ndarray, calibrations:np.ndarray, verbose:bool=True) -> tuple:

    # set calibrations
    zcal.utils.SetCalibration(scheme)
    zcal.options['verbose'] = verbose

    # Chi Squared
    def χ2(obsR:np.ndarray, obsS:np.ndarray, calR:np.ndarray, calS:np.ndarray) -> float:  
        return np.sum(np.power(obsR - calR, 2) / (obsS ** 2 + calS ** 2))

    # define prior transform
    def ptform(p:float) -> tuple:
        return (zcal.options['range'][0] - 12.) + p * ((zcal.options['range'][1] - 12.) - (zcal.options['range'][0] - 12.))
        #stats.uniform.ppf(p, loc = zcal.options['range'][0] - 12., scale = zcal.options['range'][1] - 12.)

    # define log likelihood
    def logl(uOH:float, Rs:np.ndarray, Es:np.ndarray, calibrations:np.ndarray) -> float:
        
        # verify uOH first
        if not ((uOH >= zcal.options['range'][0] - 12.) & (uOH <= zcal.options['range'][1] - 12.)):
            return -np.inf
        
        # convert to x 
        xcal = zcal.options['oh_to_x'](uOH)
        # generate all calibrated ratios
        calRs = np.array([zcal.calibrations[cal](xcal) for cal in calibrations])
        calEs = np.array([zcal.calib_errors[cal] for cal in calibrations])

        # work out and return log likelihood
        return -1 * χ2(obsR=Rs, obsS=Es, calR=calRs, calS=calEs)

    # run fitting for log(O/H)
    sampler = dynesty.NestedSampler(loglikelihood=logl, logl_args=(ratios, errors, calibrations),
                                    prior_transform=ptform,
                                    ndim=1)
    sampler.run_nested(print_progress=zcal.options['verbose'])

    # extract values from sampler
    res = sampler.results
    samples, weights = res['samples'], np.exp(res.logwt - res.logz[-1])

    # extract median and error
    ql, qm, qh = dyfunc.quantile(samples.T[0], [0.16, 0.50, 0.84], weights=weights)
    return  ql + 12., qm + 12., qh + 12.

    