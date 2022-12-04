import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel

from .conftest import test_plot_dir


def test_quad():
    params = {
        'data_dir'    : '/home/masamune/.bin/Rad1D/data',
        'tau_min'     : 1e-8,
        'tau_max'     : 1e6,
        'eps'         : 1e-4,
        'verbose'     : False,
    }

    params['n_quad'] = 2
    model = RadModel(params)
    final_iter2 = model.convergence_test()[-1]

    params['n_quad'] = 4
    model = RadModel(params)
    final_iter4 = model.convergence_test()[-1]

    params['n_quad'] = 8
    model = RadModel(params)
    final_iter8 = model.convergence_test()[-1]

    params['n_quad'] = 32
    model = RadModel(params)
    final_iter32 = model.convergence_test()[-1]

    pct_diff2 = (final_iter2 - final_iter32) * 100. / final_iter32
    pct_diff4 = (final_iter4 - final_iter32) * 100. / final_iter32
    pct_diff8 = (final_iter8 - final_iter32) * 100. / final_iter32

    fig, ax = plt.subplots(dpi=125)

    ax.plot(model.tau, pct_diff2, label='2 points')
    ax.plot(model.tau, pct_diff4, label='4 points')
    ax.plot(model.tau, pct_diff8, label='8 points')

    ax.axhline(0., ls='--', color='k')

    ax.set_xscale('log')

    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel('Percent difference')

    ax.legend()

    plt.tight_layout()
    fn = f'{test_plot_dir}/quadrature.pdf'
    fig.savefig(fn, dpi=125)
