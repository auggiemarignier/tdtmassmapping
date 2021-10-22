import numpy as np

import pytest
from utils import meanvar_from_submeanvar


@pytest.mark.flaky(reruns=10)  # rerun becuase randomness will affect accuracy
def test_meanvar_from_submeanvar():
    mu_x = 2
    var_x = 1.5
    n_x = int(1e6)
    X = np.random.normal(loc=mu_x, scale=np.sqrt(var_x), size=n_x)

    mu_y = 1
    var_y = 3
    n_y = int(7.5e6)
    Y = np.random.normal(loc=mu_y, scale=np.sqrt(var_y), size=n_y)

    Z = np.concatenate([X, Y])
    assert Z.shape[0] == n_x + n_y

    mu_z, var_z = meanvar_from_submeanvar(mu_x, mu_y, var_x, var_y, n_x, n_y)
    assert np.isclose(Z.mean(), mu_z, atol=1e-3)
    assert np.isclose(Z.std() ** 2, var_z, atol=1e-3)
