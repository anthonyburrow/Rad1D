import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel


params = {
    'data_dir'    : 'C:\dev\Rad1D\data',
    'wave_start'  : 4000.,
    'wave_end'    : 7000.,
    'n_wave'      : 500,
    'tau_max'     : 1e2,
    'eps'         : 1e-4,
    'T_eff'       : 7000,
    'n_zones'     : 256,
    'max_iter'    : 100,
    'n_quad'      : 8
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
