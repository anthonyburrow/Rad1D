from Rad1D import RadModel
import numpy as np


def test_S0_B0():
    eps_to_check = np.logspace(-4., 0.)

    for eps in eps_to_check:
        params = {
            'data_dir'    : '/home/masamune/.bin/Rad1D/data',
            'eps'         : eps,
        }

        model = RadModel(params)
        results = model.convergence_test(5000.)

        S_B = results[-1][0]
        print(f'eps = {eps} : [S(0)/B(0)] / sqrt(eps) = {S_B / np.sqrt(eps)}')

