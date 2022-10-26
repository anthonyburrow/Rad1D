import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel


params = {
    'data_dir'    : 'C:\dev\Rad1D\data',   # Absolute path to line list
    'wave_start'  : 4000.,                 # Starting wavelength
    'wave_end'    : 7000.,                 # Ending wavelength
    'cont_res'    : 0.05,                  # Points per angstrom resolved for continuum
    'line_res'    : 3.0,                   # Points per angstrom resolved for lines in line list
    'tau_max'     : 1e4,                   # Maximum tau of atmosphere for continuum
    'eps'         : 1e-4,                  # Thermalization factor
    'T_eff'       : 7000.,                 # Characteristic temperature of atmosphere
    'n_zones'     : 256,                   # Number of tau points
    'max_iter'    : 100,                   # Maximum number of lambda iterations allowed
    'eps_converge': 1e-6,                  # Factor to determine J is converged
    'n_quad'      : 8,                     # Order of Gaussian quadrature integration
}

model = RadModel(params)

synth = model.gen_spectrum(normalize=True)


# Plot spectrum
fig, ax = plt.subplots(dpi=125)

ax.plot(synth[:, 0], synth[:, 1])

ax.set_xlim(params['wave_start'], params['wave_end'])
ax.set_ylim(0., 1.05)

ax.set_xlabel('Wavelength [A]')
ax.set_ylabel('Normalized flux')

fn = './spectrum.pdf'
fig.savefig(fn, dpi=125)


# Plot T vs tau
fig, ax = plt.subplots(dpi=125)

ax.plot(model.tau, model.T)

ax.set_xlabel(r'$\tau$')
ax.set_ylabel('Temperature [K]')

fn = './T_tau.pdf'
fig.savefig(fn, dpi=125)


# Plot convergence test
lam = 5000.
results = model.convergence_test(lam)

fig, ax = plt.subplots(dpi=125)

ax.plot(model.tau[1:], results[0, 1:], c='r', ls='-', lw=1., alpha=0.4)

for i in range(1, len(results)):
    ax.plot(model.tau[1:], results[i, 1:], c='k', ls='-', lw=1.,
                     alpha=0.4)

ax.set_xscale('log')
ax.set_yscale('log')

# ax.set_xlim(-3., np.log10(params['tau_max']))
ax.set_ylim(0.05, 1.2)

ax.set_xlabel(r'$\log \tau$')
ax.set_ylabel('S / B')

plt.tight_layout()
fn = './S_B_convergence.pdf'
fig.savefig(fn, dpi=125)
