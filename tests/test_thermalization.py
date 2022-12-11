from Rad1D import RadModel
import numpy as np

from .conftest import data_dir


params = {
    'data_dir'       : data_dir,
    'Ng_accelerated' : False,
    'verbose'        : False,
}


def test_S0_B0():
    min_eps = -4.
    eps_to_check = np.logspace(min_eps, 0., int(abs(min_eps)) + 1)

    print()

    for eps in eps_to_check:
        params['max_iter'] = max(100, int(1. / eps))
        params['eps'] = eps

        model = RadModel(params)
        results = model.convergence_test()

        S_B = results[-1, 0]
        justify = int(abs(min_eps)) + 2
        print(f'  eps = {eps:<{justify}} : [S(0)/B(0)] / sqrt(eps) = {S_B / np.sqrt(eps):.6f}')
