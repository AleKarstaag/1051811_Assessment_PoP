import pytest
import numpy as np

@pytest.mark.parametrize("n,,GL_degree, A", [
    (5, 200, np.array([
        [ 4., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0.],
        [-1.,  4., -1.,  0., -1.,  0.,  0.,  0.,  0.],
        [ 0., -1.,  4.,  0.,  0., -1.,  0.,  0.,  0.],
        [-1.,  0.,  0.,  4., -1.,  0., -1.,  0.,  0.],
        [ 0., -1.,  0., -1.,  4., -1.,  0., -1.,  0.],
        [ 0.,  0., -1.,  0., -1.,  4.,  0.,  0., -1.],
        [ 0.,  0.,  0., -1.,  0.,  0.,  4., -1.,  0.],
        [ 0.,  0.,  0.,  0., -1.,  0., -1.,  4., -1.],
        [ 0.,  0.,  0.,  0.,  0., -1.,  0., -1.,  4.]])),

    (4, 300, np.array([
        [ 4., -1., -1., -0.],
        [-1.,  4., -0., -1.],
        [-1., -0.,  4., -1.],
        [-0., -1., -1.,  4.]]))
    ])

def test_bilinear_matrix(n,GL_degree, A):
    """Test insert() method."""
    from basis.grid import Poisson
    Aapprox=np.round(Poisson(n)._A(GL_degree),0)
    assert np.array_equal(Aapprox, A)

