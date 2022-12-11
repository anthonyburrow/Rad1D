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
    'Ng_accelerated' : False,
    'verbose'        : False,
}

eps_to_check = (1e-1, 1e-2, 1e-4)


def analytic_J(eps, tau):
    # In terms of B
    sqrt_3 = np.sqrt(3.)
    therm_depth = np.sqrt(3. * eps)
    J = 1. - sqrt_3 * np.exp(-therm_depth * tau) / (sqrt_3 + therm_depth)
    return J


def plot_J(ax, model, results, eps):
    if eps == 1.:
        print('J = 0 when epsilon = 1')
        return

    # Analytic J
    J_analytic = analytic_J(eps, model.tau)

    # Calculated J
    J_calc = (results[-1] - eps) / (1. - eps)

    pct_diff = 100. * (J_analytic - J_calc) / J_calc

    ax.plot(model.tau, pct_diff, '-', label=r'$\epsilon = ' + f'{eps:.0e}$')


def test_J_analytic():
    fig, ax = plt.subplots()

    for eps in eps_to_check:
        params['eps'] = eps
        model = RadModel(params)
        results = model.convergence_test(check_converged=True)

        plot_J(ax, model, results, eps)                        	

    ax.set_xscale('log')

    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel(r'Percent difference from analytic J($\tau$)')

    ax.legend()

    plt.tight_layout()
    fn = f'{test_plot_dir}/J_comparison.pdf'
    fig.savefig(fn)
