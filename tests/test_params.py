from Rad1D import RadModel
import numpy as np


def test_params():
    params = {
        'data_dir'    : '/home/masamune/.bin/Rad1D/data',
        'n_zones'     : 200,
    }

    model = RadModel(params)

    assert isinstance(model.T, np.ndarray)
    assert model.T.shape == (200, )
