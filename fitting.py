# imports
import zcal
from zcal import utils
import numpy as np
import scipy.stats as stats
import dynesty, pickle
from dynesty import utils as dyfunc
from dynesty import plotting as dyplot
from uncertainties import ufloat, unumpy as unp
from uncertainties.umath import log10 as ulog10

# main method
def FitZg_DustCorrect(id:str, scheme:str, oiii:np.ndarray, oii:np.ndarray, hb:np.ndarray, neiii:np.ndarray, correct:bool = True) -> np.ndarray:
    
    # initialise parameters
    utils.SetCalibration(scheme)

    # Generate Cardelli Grid ()
    cardelli_wl = np.linspace(0.3, 0.6, 1000)

    def logl(u:tuple) -> float:
        
        # unpack parameters
        oh, ebv = u
        
        #  convert x into correct form for calibration
        x = zcal.options['oh_to_x'](oh)
    
        if correct:
            # for each line, work out correction factors
            Aλ = utils.Cardelli_Attenuation(cardelli_wl, ebv)
            corrections = [utils.DustCorrect(zcal.wavelengths[key], wl_grid=cardelli_wl, Aλ=Aλ) for key in zcal.wavelengths]
        else:
            corrections = np.ones(len(zcal.wavelengths), dtype=float)

        # generate ufloats for objects
        if oiii[1] >= 0.:
            oiiic = ufloat((4/3)*oiii[0]*corrections[0], (4/3)*oiii[1]*corrections[0])
        else:
            oiiic = ufloat(-999.0, 0.)

        if oii[1] >= 0.:
            oiic = ufloat(oii[0]*corrections[1], oii[1]*corrections[1])
        else:
            oiic = ufloat(-999.0, 0.)

        if hb[1] >= 0.:
            hbc = ufloat(hb[0]*corrections[2], hb[1]*corrections[2])
        else:
            hbc = ufloat(-999.0, 0.)

        if neiii[1] >= 0.:
            neiiic = ufloat(neiii[0]*corrections[3], neiii[1]*corrections[3])
        else:
            neiiic = ufloat(-999.0, 0.0)

        if min(oiiic.n, oiic.n) < 0:
            return np.ones(3) * -999.0, np.ones(3) * -999.0
 
        # set up arrays
        y, yerr, calibrations = [], [], []

        # o3
        if min(oiiic.n, hbc.n) > 0:
            o3 = unp.log10(oiiic / hbc)
            y.append(unp.nominal_values(o3))
            yerr.append(unp.std_devs(o3))
            calibrations.append('O3')

        # o2
        if (min(oiic.n, hbc.n) > 0) & (scheme != 'bian'):
            o2 = unp.log10(oiic / hbc)
            y.append(unp.nominal_values(o2))
            yerr.append(unp.std_devs(o2))
            calibrations.append('O2')

        # r23
        if (min(oiiic.n, oiic.n, hbc.n) > 0) & (scheme != 'bian'):
            r23 = unp.log10((oiiic + oiic) / hbc)
            y.append(unp.nominal_values(r23))
            yerr.append(unp.std_devs(r23))
            calibrations.append('R23')

        # o32
        if min(oiiic.n, oiic.n) > 0:
            o32 = unp.log10(oiiic / oiic)
            y.append(unp.nominal_values(o32))
            yerr.append(unp.std_devs(o32))
            calibrations.append('O32')

        # ne3o2
        if min(neiiic.n, oiic.n) > 0:
            ne3o2 = unp.log10(neiiic / oiic)
            y.append(unp.nominal_values(ne3o2))
            yerr.append(unp.std_devs(ne3o2))
            calibrations.append('Ne3O2')

        # calculate model values and cast all lists as arrays
        model = np.array([zcal.calibrations[cal](x) for cal in calibrations])
        model_errs = np.array([zcal.calib_errors[cal] for cal in calibrations])
        y = np.array(y)
        yerr = np.array(yerr)

        return -0.5 * np.sum(
            (np.power(model - y, 2) / (yerr ** 2 + model_errs ** 2)) + 2 * np.log(yerr + model_errs)
        )

    def ptform(p:tuple) -> tuple:

        # unpack variables
        uOH, uEBV = p
        pOH  = (zcal.options['range'][0] - 12.) + (uOH * np.diff(zcal.options['range'])[0])

        pEBV = 2 * uEBV

        return pOH, pEBV

    sampler = dynesty.NestedSampler(loglikelihood = logl,
                                    prior_transform = ptform,
                                    ndim = 2, bootstrap = 0)

    # run sampler
    sampler.run_nested(print_progress=zcal.options['verbose'])
    results = sampler.results

    # create pickle file
    outfile = open(f'{zcal.options["res_dir"]}/samplers/{id}_{scheme}_sampler.pkl', 'wb')
    pickle.dump(sampler.results, outfile)
    outfile.close

    # extract results
    weights = np.exp(results['logwt'] - results['logz'][-1])
    oh = dyfunc.quantile(x=results['samples'][:, 0], q=[0.16, 0.50, 0.84], weights=weights)
    ebv = dyfunc.quantile(x=results['samples'][:, 1], q=[0.16, 0.50, 0.84], weights=weights)
    #print(f'12 + log(O/H) = {oh[1]+12} + {oh[2]-oh[1]} - {oh[1]-oh[0]}')
    #print(f'          EBV = {ebv[1]} + {ebv[2]-ebv[1]} - {ebv[1]-ebv[0]}')

    # make corner plot if plotting
    if zcal.options['corner_plots']:
        cfig, caxes = dyplot.cornerplot(results, color='black', labels=['12+log(O/H)', 'E(B-V)'],
                                        label_kwargs={'fontsize':25},
                                        show_titles=True)
        try:
            plot_id = int(id)
        except:
            plot_id = id
        cfig.savefig(f'{zcal.options["plot_dir"]}/corners/{plot_id}_{scheme}_corner.png')

    # return to main
    logOHp12 = np.array([oh[1]+12, oh[2]-oh[1], oh[1]-oh[0]])
    EBVs = np.array([ebv[1], ebv[2]-ebv[1], ebv[1]-ebv[0]])
    return logOHp12, EBVs

def FitZg_NoDust(id:str, scheme:str, oiii:np.ndarray, oii:np.ndarray, hb:np.ndarray, neiii:np.ndarray) -> np.ndarray:
    
    print('-> [zcal]: not implemented...')
    import sys
    sys.exit()
    # initialise parameters
    utils.SetCalibration(scheme)

    def logl(u:tuple) -> float:
        
        oh = u

        #  convert u into correct form for calibration
        x = zcal.options['oh_to_x'](oh)
    
        # generate ufloats for objects
        if oiii[1] >= 0.:
            oiiic = ufloat((4/3)*oiii[0], (4/3)*oiii[1])
        else:
            oiiic = ufloat(-999.0, 0.)

        if oii[1] >= 0.:
            oiic = ufloat(oii[0], oii[1])
        else:
            oiic = ufloat(-999.0, 0.)

        if hb[1] >= 0.:
            hbc = ufloat(hb[0], hb[1])
        else:
            hbc = ufloat(-999.0, 0.)

        if neiii[1] >= 0.:
            neiiic = ufloat(neiii[0], neiii[1])
        else:
            neiiic = ufloat(-999.0, 0.0)

        if min(oiiic.n, oiic.n) < 0:
            return np.ones(3) * -999.0, np.ones(3) * -999.0
 
        # set up arrays
        y, yerr, calibrations = [], [], []

        # o3
        if min(oiiic.n, hbc.n) > 0:
            o3 = unp.log10(oiiic / hbc)
            y.append(unp.nominal_values(o3))
            yerr.append(unp.std_devs(o3))
            calibrations.append('O3')

        # o2
        if (min(oiic.n, hbc.n) > 0) & (scheme != 'bian'):
            o2 = unp.log10(oiic / hbc)
            y.append(unp.nominal_values(o2))
            yerr.append(unp.std_devs(o2))
            calibrations.append('O2')

        # r23
        if (min(oiiic.n, oiic.n, hbc.n) > 0) & (scheme != 'bian'):
            r23 = unp.log10((oiiic + oiic) / hbc)
            y.append(unp.nominal_values(r23))
            yerr.append(unp.std_devs(r23))
            calibrations.append('R23')

        # o32
        if min(oiiic.n, oiic.n) > 0:
            o32 = unp.log10(oiiic / oiic)
            y.append(unp.nominal_values(o32))
            yerr.append(unp.std_devs(o32))
            calibrations.append('O32')

        # ne3o2
        if min(neiiic.n, oiic.n) > 0:
            ne3o2 = unp.log10(neiiic / oiic)
            y.append(unp.nominal_values(ne3o2))
            yerr.append(unp.std_devs(ne3o2))
            calibrations.append('Ne3O2')

        # calculate model values and cast all lists as arrays
        model = np.array([zcal.calibrations[cal](x) for cal in calibrations])
        model_errs = np.array([zcal.calib_errors[cal] for cal in calibrations])
        y = np.array(y)
        yerr = np.array(yerr)

        return -0.5 * np.sum(
            (np.power(model - y, 2) / (yerr ** 2 + model_errs ** 2)) + 2 * np.log(yerr + model_errs)
        )

    def ptform(p:tuple) -> tuple:
        return (zcal.options['range'][0] - 12.) + (p * np.diff(zcal.options['range'])[0])


    sampler = dynesty.NestedSampler(loglikelihood = logl,
                                    prior_transform = ptform,
                                    ndim = 1, bootstrap = 0)

    # run sampler
    sampler.run_nested(print_progress=zcal.options['verbose'])
    results = sampler.results

    # create pickle file
    outfile = open(f'./data/mosfire/{id}_{scheme}sampler.pkl', 'wb')
    pickle.dump(sampler.results, outfile)
    outfile.close

    # extract results
    weights = np.exp(results['logwt'] - results['logz'][-1])
    oh = dyfunc.quantile(x=results['samples'][:,0], q=[0.16, 0.50, 0.84], weights=weights)

    #print(f'12 + log(O/H) = {oh[1]+12} + {oh[2]-oh[1]} - {oh[1]-oh[0]}')
    #print(f'          EBV = {ebv[1]} + {ebv[2]-ebv[1]} - {ebv[1]-ebv[0]}')

    # make corner plot if plotting
    if zcal.options['corner_plots']:
        cfig, caxes = dyplot.cornerplot(results, color='black', labels=['log(O/H)',],
                                        label_kwargs={'fontsize':25},
                                        show_titles=True)
        try:
            plot_id = int(id)
        except:
            plot_id = id
        cfig.savefig(f'{zcal.options["res_dir"]}/corners/{plot_id}_{scheme}_corner_1d.png')

    # return to main
    logOHp12 = np.array([oh[1]+12, oh[2]-oh[1], oh[1]-oh[0]])
    return logOHp12