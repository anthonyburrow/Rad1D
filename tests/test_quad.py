import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel

def test_quad():
    params = {
        'data_dir'    : '/Users/adammoss/Desktop/Atmospheres/Rad1Dtest/data',   # Absolute path to line list
        'wave_start'  : 4000.,       # Starting wavelength
        'wave_end'    : 7000.,       # Ending wavelength
        'cont_res'    : 0.05,        # Points per angstrom resolved for continuum
        'line_res'    : 3.0,         # Points per angstrom resolved for lines in line list
        'tau_min'     : 1e-8,        # Minimum tau of atmosphere for continuum
        'tau_max'     : 1e6,         # Maximum tau of atmosphere for continuum
        'eps'         : 1e-4,        # Thermalization factor
        'T_eff'       : 6000.,       # Characteristic temperature of atmosphere
        'n_zones'     : 256,         # Number of tau points
        'max_iter'    : 200,         # Maximum number of lambda iterations allowed
        'accelerated' : True,        # Use ALI algorithm for faster convergence
        'eps_converge': 1e-7,        # Factor to determine J is converged
        'n_quad'      : 32,          # Order of Gaussian quadrature integration
        'verbose'     : True,        # Display stdout output
    }

    params['n_quad'] = 2
    model = RadModel(params)
    two_results = model.convergence_test()
    two_final = two_results[-1]

    params['n_quad'] = 4
    model = RadModel(params)
    four_results = model.convergence_test()
    four_final = four_results[-1]

    params['n_quad'] = 8
    model = RadModel(params)
    eight_results = model.convergence_test()
    eight_final = eight_results[-1]

    params['n_quad'] = 32
    model = RadModel(params)
    thirtytwo_results = model.convergence_test()
    thirtytwo_final = thirtytwo_results[-1]

    two_res=[]
    four_res=[]
    eight_res=[]
    two_perc=[]
    four_perc=[]
    eight_perc=[]
    for i in range(len(thirtytwo_final)):
        two_res.append(two_final[i]-thirtytwo_final[i])
        four_res.append(four_final[i]-thirtytwo_final[i])
        eight_res.append(eight_final[i]-thirtytwo_final[i])
        two_perc.append((two_res[i]/thirtytwo_final[i])*100)
        four_perc.append((four_res[i]/thirtytwo_final[i])*100)
        eight_perc.append((eight_res[i]/thirtytwo_final[i])*100)

    fig, ax = plt.subplots(dpi=125)
    ax.plot(model.tau, two_final, ls='-', lw=1., alpha=0.4,label='2')
    ax.plot(model.tau, four_final, ls='-', lw=1., alpha=0.4,label='4')
    ax.plot(model.tau, eight_final, ls='-', lw=1., alpha=0.4,label='8')
    ax.plot(model.tau, thirtytwo_final, ls='-', lw=1., alpha=0.4,label='32')

    ax.set_xscale('log')
    ax.set_yscale('log')

    #ax.set_xlim(-3., np.log10(params['tau_max']))
    #ax.set_ylim(0.9, 1.2)

    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel('S / B')
    ax.legend(loc='best')
    plt.tight_layout()
    fn = './quad.png'
    fig.savefig(fn, dpi=125)
    
    fig, ax = plt.subplots(dpi=125)
    ax.plot(model.tau, two_res, label='2')
    ax.plot(model.tau, four_res, label='4')
    ax.plot(model.tau, eight_res, label='8')
    ax.set_xscale('log')
    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel('Residual')
    ax.legend(loc='best')
    plt.tight_layout()
    fn = './resid.png'
    fig.savefig(fn, dpi=125)
    
    fig, ax = plt.subplots(dpi=125)
    ax.plot(model.tau, two_perc, label='2')
    ax.plot(model.tau, four_perc, label='4')
    ax.plot(model.tau, eight_perc, label='8')
    ax.set_xscale('log')
    ax.set_xlabel(r'$\tau$')
    ax.set_ylabel('% Difference')
    ax.legend(loc='best')
    ax.axhline(0, ls='--', color='r')
    plt.tight_layout()
    fn = './percdiff.png'
    fig.savefig(fn, dpi=125)
