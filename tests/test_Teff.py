from Rad1D import RadModel
import numpy as np
import matplotlib.pyplot as plt

from .conftest import test_plot_dir, data_dir
from .plot_setup import paper_plot
paper_plot()


params = {
    'data_dir'       : data_dir,
    'wave_start'     : 3000.,
    'wave_end'       : 7000.,
    'T_eff'          : 6000.,
    'Ng_accelerated' : False,
    'verbose'        : False,
}


def get_T_wien(lam):
    return 2.8977719e7 / lam


def test_Teff():
    print()

    model = RadModel(params)
    spectrum = model.gen_spectrum()

    print(f'  Input temperature: {params["T_eff"]:.3f} K')

    max_ind = spectrum[:, 1].argmax()
    model_Teff = get_T_wien(spectrum[max_ind, 0])
    tau_ind = np.abs(model.T - model_Teff).argmin()
    closest_tau = model.tau[tau_ind]

    print(f'  Spectrum BB temperature: {model_Teff:.3f} K (tau = ~{closest_tau:.4f})')

    tau_eff = 0.64
    expected_ind = np.abs(model.tau - tau_eff).argmin()
    expected_Teff = model.T[expected_ind]

    print(f'  Expected temperature (tau = {tau_eff:.2f}): {expected_Teff:.3f} K')

    # Plot spectra
    fig, ax = plt.subplots(dpi=125)

    ax.plot(spectrum[:,0], spectrum[:,1], label='input')

    ax.set_xlim(params['wave_start'], params['wave_end'])

    ax.set_xlabel(r'Wavelength [$\AA$]')
    ax.set_ylabel('Emergent flux')

    plt.tight_layout()
    fn = f'{test_plot_dir}/test_spectrum.pdf'
    fig.savefig(fn, dpi=125)
