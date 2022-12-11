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


params['n_quad'] = 32
model = RadModel(params)
final_iter32 = model.convergence_test()[-1]


def plot_pct_diff(ax, quad):
    params['n_quad'] = quad
    model = RadModel(params)
    final_iter = model.convergence_test()[-1]

    pct_diff = 100. * (final_iter - final_iter32) / final_iter32

    ax.plot(model.tau, pct_diff, label=fr'$S_{quad}$')


def test_quad():
    fig, ax = plt.subplots(dpi=125)

    plot_pct_diff(ax, 2)
    plot_pct_diff(ax, 4)
    plot_pct_diff(ax, 8)

    ax.axhline(0., ls='--', color='k')

    ax.set_xscale('log')

    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel(r'% Difference from $S_{32}$')

    ax.legend()

    plt.tight_layout()
    fn = f'{test_plot_dir}/quadrature.pdf'
    fig.savefig(fn, dpi=125)
