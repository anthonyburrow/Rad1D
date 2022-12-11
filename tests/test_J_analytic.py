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

hc = 12398.4198
c = 29979.2458
k = 8.617333262e-5
k_hc = k / hc


def planck(lam, T):
    bbScale = 2. * hc * c * 1.e8
    bb = bbScale / (pow(lam, 5.) * (np.exp(1. / (k_hc * T * lam)) - 1.))
    return bb


def analytic_J(eps, tau, a=0., b=0.):
    c = np.sqrt(3. * eps)
    J = a + (b - c * a) * np.exp(-c * tau) / (np.sqrt(3.) + c)
    return J


def plot_J(ax, model, results, eps):
    if eps == 1.:
        return

    # Find analytic J
    a = planck(5000., params['T_eff'])
    J_analytic = analytic_J(eps, model.tau, a=a)

    # Find calculated J
    J_calc = (results - eps * a) / (1. - eps)
    diffs = np.abs(J_calc[1:] - J_calc[:-1]) / J_calc[:-1]
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

    J_calc = J_calc[ind_needed]

    # Plot here
    pct_diff = 100. * (J_analytic - J_calc) / J_calc

    ax.plot(model.tau, pct_diff, '-', label=r'$\epsilon = ' + f'{eps:.0e}$')


def test_J_analytic():
    fig, ax = plt.subplots()

    for i, eps in enumerate(eps_to_check):
        params['eps'] = eps
        model = RadModel(params)
        results = model.convergence_test()

        plot_J(ax, model, results, eps)                        	

    ax.set_xscale('log')

    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel('% Difference')

    ax.legend()

    plt.tight_layout()
    fn = f'{test_plot_dir}/J_comparison.pdf'
    fig.savefig(fn)
