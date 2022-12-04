import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel

from .conftest import test_plot_dir


def test_epsilon():
    params = {
        'data_dir'    : '/Users/adammoss/Desktop/Atmospheres/Rad1Dtest/data',
        'wave_start'  : 4000.,
        'wave_end'    : 7000.,
        'cont_res'    : 0.05,
        'line_res'    : 3.0,
        'tau_min'     : 1e-8,
        'tau_max'     : 1e6,
        'eps'         : 1e-4,
        'T_eff'       : 6000.,
        'n_zones'     : 256,
        'max_iter'    : 1000,
        'accelerated' : True,
        'eps_converge': 1e-7,
        'n_quad'      : 32,
        'verbose'     : True,
    }

    n_iter=params['max_iter']
    n_zones=params['n_zones']
    params['eps'] = 1e-2
    eps=params['eps']
    eps_converge=params['eps_converge']
    eps_con_two = np.sqrt(1e-2)
    model = RadModel(params)
    two_results = model.convergence_test()
    J = (two_results - eps) / (1 - eps)
    diff = np.abs(J[1:] - J[:-1]) / J[:-1]
    k=0
    for i in range(n_iter-1):
        for j in range(n_zones):
            if diff[i,j] < eps_converge:
               k+=1
            if k == n_zones:
                print(i)
                two_correct = two_results[0:i]
                break
        if k == n_zones:
            break
        k=0 

    params['eps'] = 1e-4
    eps=params['eps']
    eps_con_four = np.sqrt(1e-4)
    model = RadModel(params)
    four_results = model.convergence_test()
    J = (four_results - eps) / (1 - eps)
    diff = np.abs(J[1:] - J[:-1]) / J[:-1]
    k=0
    for i in range(n_iter-1):
        for j in range(n_zones):
            if diff[i,j] < eps_converge:
               k+=1
            if k == n_zones:
                print(i)
                four_correct = four_results[0:i]
                break
        if k == n_zones:
            break
        k=0  

    params['eps'] = 1e-6
    eps=params['eps']
    eps_con_six = np.sqrt(1e-6)
    model = RadModel(params)
    six_results = model.convergence_test()
    J = (six_results - eps) / (1 - eps)
    diff = np.abs(J[1:] - J[:-1]) / J[:-1]
    k=0
    for i in range(n_iter-1):
        for j in range(n_zones):
            if diff[i,j] < eps_converge:
               k+=1
            if k == n_zones:
                print(i)
                six_correct = six_results[0:i]
                break
        if k == n_zones:
            break
        k=0 


    fig, ax = plt.subplots(dpi=125)
    ax.plot(model.tau[1:], two_correct[0, 1:], c='r', ls='-', lw=1., alpha=0.4,label='1e-2,115 iterations')
    for i in range(1, len(two_correct)):
        ax.plot(model.tau[1:], two_correct[i, 1:], c='k', ls='-', lw=1.,
                     alpha=0.4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel('S/B')
    ax.legend(loc='best')
    ax.axhline(eps_con_two, ls='--', color='r')
    plt.tight_layout()
    fn = f'{test_plot_dir}/eps2.png'
    fig.savefig(fn, dpi=125)
  
    fig, ax = plt.subplots(dpi=125)
    ax.plot(model.tau[1:], four_correct[0, 1:], c='r', ls='-', lw=1., alpha=0.4, label='1e-4, 372 iterations')
    for i in range(1, len(four_correct)):
        ax.plot(model.tau[1:], four_correct[i, 1:], c='k', ls='-', lw=1.,
                     alpha=0.4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel('S/B')
    ax.legend(loc='best')
    ax.axhline(eps_con_four, ls='--', color='r')
    plt.tight_layout()
    fn = f'{test_plot_dir}/eps4.png'
    fig.savefig(fn, dpi=125)
  
    fig, ax = plt.subplots(dpi=125)
    ax.plot(model.tau[1:], six_correct[0, 1:], c='r', ls='-', lw=1., alpha=0.4,label='1e-6, 605 iterations')
    for i in range(1, len(six_correct)):
        ax.plot(model.tau[1:], six_correct[i, 1:], c='k', ls='-', lw=1.,
                     alpha=0.4)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel('S/B')
    ax.legend(loc='best')
    ax.axhline(eps_con_six, ls='--', color='r')
    plt.tight_layout()
    fn = f'{test_plot_dir}/eps6.png'
    fig.savefig(fn, dpi=125)
                              	
