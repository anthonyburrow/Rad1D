from Rad1D import RadModel
import numpy as np


def get_T_wien(lam):
    return 2.8977719e7 / lam


def test_params():
    print()

    params = {
        'data_dir'    : '/home/masamune/.bin/Rad1D/data',
        'T_eff'       : 7000.,
        'verbose'     : False,
    }

    model = RadModel(params)

    spectrum = model.gen_spectrum(normalize=True)

    print(f'  Input temperature: {params["T_eff"]:.3f} K')

    max_ind = spectrum[:, 1].argmax()
    model_Teff = get_T_wien(spectrum[max_ind, 0])
    print(f'  Spectrum BB temperature: {model_Teff:.3f} K')

    tau_eff = 0.64
    expected_ind = abs(model.tau - tau_eff).argmin()
    expected_Teff = model.T[expected_ind]
    print(f'  Expected temperature (tau = {tau_eff:.2f}): {expected_Teff:.3f} K')
