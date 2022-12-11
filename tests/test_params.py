from Rad1D import RadModel
import numpy as np

from .conftest import data_dir


def test_params():
    params = {
        'data_dir'    : data_dir,
        'n_zones'     : 200,
        'verbose'     : False,
    }

    model = RadModel(params)

    assert isinstance(model.T, np.ndarray)
    assert model.T.shape == (200, )
