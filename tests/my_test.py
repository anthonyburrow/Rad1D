import numpy as np
import matplotlib.pyplot as plt
from Rad1D import RadModel


params = {
    'data_dir'    : 'C:\dev\Rad1D\data',
    'out_filename': 'synthetic.dat',
    'wave_start'  : 4000.,
    'wave_end'    : 7000.,
    'n_wave'      : 1000,
    'tau_max'     : 1000,
    'eps'         : 0.001,
    'T_eff'       : 6000,
    'n_zones'     : 256,
    'max_iter'    : 100,
    'n_quad'      : 8
}

model = RadModel(params)

# synth = model.gen_spectrum()
# synth[:, 1] /= synth[:, 1].max()

# Plot spectrum
# fig, ax = plt.subplots(dpi=125)

# ax.plot(synth[:, 0], synth[:, 1])

# ax.set_xlabel('Wavelength [A]')
# ax.set_ylabel('Normalized flux')

# fn = './spectrum.pdf'
# fig.savefig(fn, dpi=125)


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
print(results)

fig, ax = plt.subplots(dpi=125)

for i in range(len(results)):
    ax.plot(np.log10(model.tau[1:]), results[i, 1:], c='k', ls='-', lw=1.)

ax.set_xlim(-3., np.log10(params['tau_max']))
ax.set_ylim(0., 1.05)

ax.set_xlabel(r'$\log \tau$')
ax.set_ylabel('S / B')

fn = './S_B_convergence.pdf'
fig.savefig(fn, dpi=125)
