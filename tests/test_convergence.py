import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel

from .conftest import test_plot_dir


def test_convergence():
    params = {
        'data_dir'    : '/home/masamune/.bin/Rad1D/data',
        'eps'         : 1e-4,
        'eps_converge' : 1e-8,
        'max_iter'    : 1000,
        'verbose'     : False,
    }
    eps = params['eps']
    eps_converge = params['eps_converge']
    max_iter = params['max_iter']
    model = RadModel(params)
    Ng_surface = model.convergence_test()[:, 0]
    Ng_results = model.convergence_test()    

    J_Ng = (Ng_results - eps) / (1 - eps)
    diffs_Ng = np.abs(J_Ng[1:] - J_Ng[:-1]) / J_Ng[:-1]
    max_diffs_Ng = diffs_Ng.max(axis=1)

#    for i in range(max_iter):
#        if max_diffs_Ng[i] > eps_converge:
#            continue
#        ind_needed = i
#        break
#    else:
#        ind_needed = max_iter
#    Ng_results = Ng_results[:ind_needed]

    params['Ng_accelerated'] = False
    model = RadModel(params)
    ali_surface = model.convergence_test()[:, 0]
    ali_results = model.convergence_test()


    J_ali = (ali_results - eps) / (1 - eps)
    diffs_ali = np.abs(J_ali[1:] - J_ali[:-1]) / J_ali[:-1]
    max_diffs_ali = diffs_ali.max(axis=1)

#    for i in range(max_iter):
#        if max_diffs_ali[i] > eps_converge:
#            continue
#        ind_needed = i
#        break
#    else:
#        ind_needed = max_iter
#    ali_results = ali_results[:ind_needed]    

    params['accelerated'] = False
    model = RadModel(params)
    li_surface = model.convergence_test()[:, 0]
    li_results = model.convergence_test()

    J_li = (li_results - eps) / (1 - eps)
    diffs_li = np.abs(J_li[1:] - J_li[:-1]) / J_li[:-1]
    max_diffs_li = diffs_li.max(axis=1)

#    for i in range(max_iter):
#        if max_diffs_li[i] > eps_converge:
#            continue
#        ind_needed = i
#        break
#    else:
#        ind_needed = max_iter
#    li_results = li_results[:ind_needed]

    iterations = np.arange(1, params['max_iter'] + 1)
    iterations_short = np.arange(1, params['max_iter'])
    eps = np.sqrt(params['eps'])

    fig3, ax3 = plt.subplots(1, 2, figsize=(12, 4), sharex=True)
    ax3[0].plot(iterations[:len(Ng_surface)], Ng_surface, color='tab:orange', label='Ng')
    ax3[0].plot(iterations, ali_surface, color='tab:blue', label='ALI')
    ax3[0].plot(iterations, li_surface, color='k', label='LI')

    ax3[0].axhline(eps, ls='--')

    ax3[0].set_xlabel('Iteration')
    ax3[0].set_ylabel(r'$(S / B)[0]$')

    ax3[0].set_xscale('log')
    ax3[0].set_yscale('log')

    ax3[0].legend(loc='best')
   
    ax3[1].plot(iterations_short, max_diffs_Ng, color='tab:orange', label='Ng')
    ax3[1].plot(iterations_short, max_diffs_ali, color='tab:blue', label='ALI')
    ax3[1].plot(iterations_short, max_diffs_li, color='k', label='LI')


    ax3[1].set_xlabel('Iteration')
    ax3[1].set_ylabel(r'$\Delta J/J$')

    ax3[1].set_xscale('log')
    ax3[1].set_yscale('log')

#    ax3[1].legend(loc='best')

    plt.tight_layout()

    fn3 = f'{test_plot_dir}/iterations.pdf'
    fig3.savefig(fn3, dpi=125)

