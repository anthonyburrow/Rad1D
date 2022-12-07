import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel

from .conftest import test_plot_dir


def test_convergence():
    params = {
        'data_dir'    : '/home/masamune/.bin/Rad1D/data',
        'eps'         : 1e-4,
        'max_iter'    : 1000,
        'verbose'     : False,
    }
    model = RadModel(params)
    Ng_result = model.convergence_test()[:, 0]

    params['Ng_accelerated'] = False
    model = RadModel(params)
    ali_result = model.convergence_test()[:, 0]

    params['accelerated'] = False
    model = RadModel(params)
    li_result = model.convergence_test()[:, 0]

    iterations = np.arange(1, params['max_iter'] + 1)

    fig, ax = plt.subplots(dpi=125)

    ax.plot(iterations[:len(Ng_result)], Ng_result, color='tab:orange', label='Ng')
    ax.plot(iterations, ali_result, color='tab:blue', label='ALI')
    ax.plot(iterations, li_result, color='k', label='LI')

    eps = np.sqrt(params['eps'])
    ax.axhline(eps, ls='--')

    ax.set_xlabel('Iteration')
    ax.set_ylabel(r'$(S / B)[0]$')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend()

    plt.tight_layout()

    fn = f'{test_plot_dir}/iterations.pdf'
    fig.savefig(fn, dpi=125)
