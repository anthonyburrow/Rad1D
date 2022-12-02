import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel


def test_convergence():
    params = {
        'data_dir'    : '/home/masamune/.bin/Rad1D/data',
        'eps'         : 1e-6,
        'max_iter'    : 1000,
        'verbose'     : False,
    }

    params['accelerated'] = True
    model = RadModel(params)
    ali_result = model.convergence_test()[:, 0]

    params['accelerated'] = False
    model = RadModel(params)
    li_result = model.convergence_test()[:, 0]

    iterations = np.arange(1, params['max_iter'] + 1)

    fig, ax = plt.subplots(dpi=125)

    ax.plot(iterations, ali_result, color='r', label='ALI')
    ax.plot(iterations, li_result, color='k', label='LI')

    eps = np.sqrt(params['eps'])
    ax.axhline(eps, ls='--')

    ax.set_xlabel('Iteration')
    ax.set_ylabel(r'$(S / B)[0]$')
    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend()

    plt.tight_layout()

    fn = './iterations.pdf'
    fig.savefig(fn, dpi=125)