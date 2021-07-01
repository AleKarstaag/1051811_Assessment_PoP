import pytest
import numpy as np
import numpy.testing as npt

@pytest.mark.parametrize("n, A", [
    (5, np.array([
        [ 4., -1.,  0., -1.,  0.,  0.,  0.,  0.,  0.],
        [-1.,  4., -1.,  0., -1.,  0.,  0.,  0.,  0.],
        [ 0., -1.,  4.,  0.,  0., -1.,  0.,  0.,  0.],
        [-1.,  0.,  0.,  4., -1.,  0., -1.,  0.,  0.],
        [ 0., -1.,  0., -1.,  4., -1.,  0., -1.,  0.],
        [ 0.,  0., -1.,  0., -1.,  4.,  0.,  0., -1.],
        [ 0.,  0.,  0., -1.,  0.,  0.,  4., -1.,  0.],
        [ 0.,  0.,  0.,  0., -1.,  0., -1.,  4., -1.],
        [ 0.,  0.,  0.,  0.,  0., -1.,  0., -1.,  4.]])),

    (4, np.array([
        [ 4., -1., -1., -0.],
        [-1.,  4., -0., -1.],
        [-1., -0.,  4., -1.],
        [-0., -1., -1.,  4.]])),
    
    (6, np.array([
        [ 4., -1.,  0.,  0., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.],
        [-1.,  4., -1.,  0.,  0., -1.,  0.,  0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.],
        [ 0., -1.,  4., -1.,  0.,  0., -1., -0.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.],
        [ 0.,  0., -1.,  4.,  0.,  0., -0., -1.,  0.,  0.,  0.,  0.,  0., 0.,  0.,  0.],
        [-1.,  0.,  0.,  0.,  4., -1.,  0.,  0., -1.,  0.,  0.,  0.,  0., 0.,  0.,  0.],
        [ 0., -1.,  0.,  0., -1.,  4., -1.,  0.,  0., -1., -0.,  0.,  0., 0.,  0.,  0.],
        [ 0.,  0., -1.,  0.,  0., -1.,  4., -1.,  0., -0., -1.,  0.,  0., 0.,  0.,  0.],
        [ 0.,  0., -0., -1.,  0.,  0., -1.,  4.,  0.,  0.,  0., -1.,  0., 0.,  0.,  0.],
        [ 0.,  0.,  0.,  0., -1.,  0.,  0.,  0.,  4., -1.,  0.,  0., -1., 0.,  0.,  0.],
        [ 0.,  0.,  0.,  0.,  0., -1., -0.,  0., -1.,  4., -1.,  0.,  0.,-1.,  0.,  0.],
        [ 0.,  0.,  0.,  0.,  0., -0., -1.,  0.,  0., -1.,  4., -1.,  0., 0., -1.,  0.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  0.,  0., -1.,  4.,  0., 0.,  0., -1.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1., -0.,  0.,  0.,  4.,-1.,  0.,  0.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -0., -1.,  0.,  0., -1., 4., -1.,  0.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  0.,  0.,-1.,  4., -1.],
        [ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0., -1.,  0., 0., -1.,  4.]]))
    ])

def test_stiffness_matrix(n, A):

    from basis import Poisson
    P=Poisson(n)
    P._A()
    Aapprox=P.A
    npt.assert_almost_equal(Aapprox, A, decimal=0)

@pytest.mark.parametrize("f,length,n,boundary_value",[
    (lambda x,y: np.sin(x+y),100,5,np.zeros(100)),

    (lambda x,y: x**2+y**2,200,4,np.zeros(200)),

    (lambda x,y: np.cos(x-y)+np.log(x),300,6,np.zeros(300))
])

def test_boundary_condition(f,length,n,boundary_value):
    from basis import Poisson
    P=Poisson(n,length,[0,0],f)
    P.uh(50,50)
    x0=boundary_value
    y0=boundary_value
    xlength=boundary_value
    ylength=boundary_value

    for i in range(length):
        
        x0[i]=np.round(P.uh(0,i),17)
        y0[i]=np.round(P.uh(i,0),17)
        xlength[i]=np.round(P.uh(length,i),17)
        ylength[i]=np.round(P.uh(i,length),17)
    
    assert np.array_equal(x0,boundary_value)
    assert np.array_equal(y0,boundary_value)
    assert np.array_equal(xlength,boundary_value)
    assert np.array_equal(ylength,boundary_value)


       
    