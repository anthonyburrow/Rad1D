import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel

from .conftest import test_plot_dir, data_dir
from .plot_setup import paper_plot
paper_plot()


params = {
    'data_dir'       : data_dir,
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
    iter_needed = len(results)

    if eps == 1.:
        ax.plot(model.tau, results[0], c='r', ls='-', lw=1., alpha=0.4)
        for i in range(1, iter_needed):
            ax.plot(model.tau, results[i], c='k', ls='-', lw=1., alpha=0.4)
        ax.set_title(r'$\epsilon = 1.0$')
        return

    ax.plot(model.tau, results[0], c='r', ls='-', lw=1., alpha=0.4)
    for i in range(1, iter_needed):
        ax.plot(model.tau, results[i], c='k', ls='-', lw=1., alpha=0.4)

    ax.axhline(np.sqrt(eps), ls='--', color='k')

    ax.text(0.60, 0.1, f'{iter_needed} iterations', transform=ax.transAxes)
    ax.set_title(r'$\epsilon = $' + f'{eps:.0e}')


def test_epsilon():
    fig, ax = plt.subplots(1, 3, figsize=(11, 4), sharey=True)

    for i, eps in enumerate(eps_to_check):
        params['eps'] = eps
        model = RadModel(params)
        results = model.convergence_test(check_converged=True)

        plot_eps(ax[i], model, results, eps)

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
