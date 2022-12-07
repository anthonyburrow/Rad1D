from Rad1D import RadModel
import numpy as np
import matplotlib.pyplot as plt

from .conftest import test_plot_dir


def get_T_wien(lam):
    return 2.8977719e7 / lam


def test_Teff():
    print()

    params = {
        'data_dir'    : '/home/masamune/.bin/Rad1D/data',
        'T_eff'       : 6000.,
        'verbose'     : False,
        'wave_start'  : 3000.,
        'wave_end'    : 7000.,
    }
   
    model = RadModel(params)

    spectrum = model.gen_spectrum(normalize=False)
  
    print(f'  Input temperature: {params["T_eff"]:.3f} K')

    max_ind = spectrum[:, 1].argmax()
    model_Teff = get_T_wien(spectrum[max_ind, 0])
    tau_ind = np.abs(model.T - model_Teff).argmin()
    closest_tau = model.tau[tau_ind]
    print(f'  Spectrum BB temperature: {model_Teff:.3f} K (tau = ~{closest_tau:.4f})')
    
    params['T_eff']=model_Teff
    BBmodel = RadModel(params)
    BBspectrum = BBmodel.gen_spectrum(normalize=False)

    tau_eff = 0.64
    expected_ind = np.abs(model.tau - tau_eff).argmin()
    expected_Teff = model.T[expected_ind]
    print(f'  Expected temperature (tau = {tau_eff:.2f}): {expected_Teff:.3f} K')

    params['T_eff']=expected_Teff
    expectmodel = RadModel(params)
    expectspectrum = expectmodel.gen_spectrum(normalize=False)

    # Plot spectra
    fig, ax = plt.subplots(dpi=125)

    ax.plot(spectrum[:,0],spectrum[:,1],label='input')
    ax.plot(BBspectrum[:,0],BBspectrum[:,1],label='Wien')
#    ax.plot(expectspectrum[:,0],expectspectrum[:,1],label='tau=2/3')

    ax.set_xlim(params['wave_start'], params['wave_end'])
#    ax.set_ylim(0., 1.05)

    ax.set_xlabel('Wavelength [A]')
    ax.set_ylabel('Flux')
    ax.legend(loc='best')

    fn = f'{test_plot_dir}/test_spectrum.pdf'
    fig.savefig(fn, dpi=125)
