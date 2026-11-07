def test_package_imports():
    import rovsim.sim
    assert rovsim.sim is not None


def test_numpy_available():
    import numpy as np
    assert int(np.__version__.split(".")[0]) >= 2
