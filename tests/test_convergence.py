import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel

from .conftest import test_plot_dir


params = {
    'data_dir'     : '/home/masamune/.bin/Rad1D/data',
    'eps'          : 1e-4,
    'eps_converge' : 1e-8,
    'max_iter'     : 1000,
    'verbose'      : False,
}

iterations = np.arange(1, params['max_iter'] + 1)


def plot_SoverB(ax, results, *args, **kwargs):
    surface_results = results[:, 0]

    # Plot
    ax.plot(iterations, surface_results, *args, **kwargs)


def plot_deltaJ(ax, results, *args, **kwargs):
    eps = params['eps']

    J = (results - eps) / (1. - eps)
    diffs = 2. * np.abs(J[1:] - J[:-1]) / (J[1:] + J[:-1])
    max_diffs = diffs.max(axis=1)

    ax.plot(iterations[1:], max_diffs, *args, **kwargs)


def test_convergence():
    eps = params['eps']
    eps_converge = params['eps_converge']

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))

    # Ng accleration
    model = RadModel(params)
    results = model.convergence_test(check_converged=False)

    plot_SoverB(ax[0], results, label='Ng')
    plot_deltaJ(ax[1], results)

    # ALI
    params['Ng_accelerated'] = False
    model = RadModel(params)
    results = model.convergence_test(check_converged=False)

    plot_SoverB(ax[0], results, label='ALI')
    plot_deltaJ(ax[1], results)

    # Lambda iteration
    params['accelerated'] = False
    model = RadModel(params)
    results = model.convergence_test(check_converged=False)

    plot_SoverB(ax[0], results, label='LI')
    plot_deltaJ(ax[1], results)

    # Plot
    ax[0].axhline(np.sqrt(eps), ls='--', c='k')
    ax[1].axhline(eps_converge, ls='--', c='k')

    [_ax.set_xscale('log') for _ax in ax]
    [_ax.set_yscale('log') for _ax in ax]

    [_ax.set_xlabel('Iteration') for _ax in ax]
    ax[0].set_ylabel(r'$(S / B)[0]$')
    ax[1].set_ylabel(r'max($\Delta J$ / $J$)')

    [_ax.set_xlim(1, params['max_iter']) for _ax in ax]
    ax[0].set_ylim(top=1.)
    ax[1].set_ylim(top=1.)

    ax[0].legend()  

    plt.tight_layout()

    fn = f'{test_plot_dir}/iterations.pdf'
    fig.savefig(fn, dpi=125)

