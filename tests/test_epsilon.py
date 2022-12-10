import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel

from .conftest import test_plot_dir
from .plot_setup import paper_plot
paper_plot()


params = {
    'data_dir'       : '/Users/adammoss/Desktop/Atmospheres/Rad1Dtest/data',
    'tau_min'        : 1e-8,
    'tau_max'        : 1e6,
    'T_eff'          : 6000.,
    'max_iter'       : 1000,
    'eps_converge'   : 1e-8,
    'Ng_accelerated' : False,
    'verbose'        : False,
}

eps_to_check = (1e0, 1e-1, 1e-4)


def plot_eps(ax, model, results, eps):
    if eps == 1.:
        ax.plot(model.tau, results[0], c='r', ls='-', lw=1., alpha=0.4)
        for i in range(1, params['max_iter']):
            ax.plot(model.tau, results[i], c='k', ls='-', lw=1., alpha=0.4)
        ax.set_title(r'$\epsilon = 1.0$')
        return

    J = (results - eps) / (1 - eps)
    diffs = 2. * np.abs(J[1:] - J[:-1]) / (J[1:] + J[:-1])
    max_diffs = diffs.max(axis=1)

    eps_converge = params['eps_converge']
    max_iter = params['max_iter']
    for i in range(max_iter):
        if max_diffs[i] > eps_converge:
            continue
        ind_needed = i
        break
    else:
        ind_needed = max_iter

    results = results[:ind_needed]

    ax.plot(model.tau, results[0], c='r', ls='-', lw=1., alpha=0.4)
    for i in range(1, ind_needed):
        ax.plot(model.tau, results[i], c='k', ls='-', lw=1., alpha=0.4)

    ax.axhline(np.sqrt(eps), ls='--', color='k')

    ax.text(0.60, 0.1, f'{ind_needed + 1} iterations', transform=ax.transAxes)
    ax.set_title(r'$\epsilon = $' + f'{eps:.0e}')


def test_epsilon():
    fig, ax = plt.subplots(1, 3, figsize=(11, 4), sharey=True)

    for i, eps in enumerate(eps_to_check):
        params['eps'] = eps
        model = RadModel(params)
        results = model.convergence_test()

        plot_eps(ax[i], model, results, eps)

    # Fig 1 properties
    [_ax.set_xlim(params['tau_min'], params['tau_max']) for _ax in ax]
    bottom_convergence = np.sqrt(np.array(eps_to_check).min())
    ax[0].set_ylim(0.75 * bottom_convergence, 1.25)

    [_ax.set_xscale('log') for _ax in ax]
    ax[0].set_yscale('log')

    [_ax.set_xlabel(r'$\tau$') for _ax in ax]
    ax[0].set_ylabel('S / B')

    plt.tight_layout()
    fn = f'{test_plot_dir}/eps_convergence.pdf'
    fig.savefig(fn)                              	
