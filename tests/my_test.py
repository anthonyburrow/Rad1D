import matplotlib.pyplot as plt
from Rad1D import RadModel


params = {
    'data_dir'    : 'C:\dev\Rad1D\data',
    'out_filename': 'synthetic.dat',
    'wave_start'  : 4000.,
    'wave_end'    : 7000.,
    'n_wave'      : 1000,
    'tau_max'     : 100,
    'eps'         : 0.1,
    'T_eff'         : 6000,
    'n_zones'     : 256,
    'max_iter'    : 100,
    'n_quad'      : 8
}

model = RadModel(params)

synth = model.gen_spectrum()
synth[:, 1] /= synth[:, 1].max()

# Plot spectrum
fig, ax = plt.subplots(dpi=125)

ax.plot(synth[:, 0], synth[:, 1])

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
