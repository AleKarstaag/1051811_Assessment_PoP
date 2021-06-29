import pytest
import numpy as np
import numpy.testing as npt

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
        [-0., -1., -1.,  4.]])),
    
    (6,100,np.array([
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

def test_bilinear_matrix(n,GL_degree, A):
    """Test insert() method."""
    from basis import Poisson
    Aapprox=Poisson(n)._A(GL_degree)
    npt.assert_almost_equal(Aapprox, A, decimal=0)

@pytest.mark.parametrize("f,length,n,GLdegree,boundary_value",[
    (lambda x,y: np.sin(x+y),100,5,200,np.zeros(100)),

    (lambda x,y: x**2+y**2,200,4,300,np.zeros(200)),

    (lambda x,y: np.cos(x-y)+np.log(x),300,6,100,np.zeros(300))
])

def test_closure(f,length,n,GLdegree,boundary_value):
    from basis import Poisson
    P=Poisson(n,length)
    P._A(GLdegree)
    P._L(f,GLdegree)
    P._U()

    x0=boundary_value
    y0=boundary_value
    xlength=boundary_value
    ylength=boundary_value

    for i in range(length):
        x0[i]=np.round(P.u(0,i),1)
        y0[i]=np.round(P.u(i,0),1)
        xlength[i]=np.round(P.u(length,i),1)
        ylength[i]=np.round(P.u(i,length),1)
    
    assert np.array_equal(x0,boundary_value)
    assert np.array_equal(y0,boundary_value)
    assert np.array_equal(xlength,boundary_value)
    assert np.array_equal(ylength,boundary_value)


@pytest.mark.parametrize("f,length,n,GLdegree",[
    (lambda x,y: np.sin(x+y),100,5,200),

    (lambda x,y: x**2+y**2,200,4,300),

    (lambda x,y: np.cos(x-y)+np.log(x),300,6,100)
])

def test_AttributeError(f, length, n, GLdegree):
    from basis import Poisson
    P=Poisson(n,length)
    with pytest.raises(AttributeError):
        P._U()

@pytest.mark.parametrize("f,length,n,GLdegree",[
    (lambda x,y: np.sin(x+y),100,5,200),

    (lambda x,y: x**2+y**2,200,4,300),

    (lambda x,y: np.cos(x-y)+np.log(x),300,6,100)
])

def test_AttributeError2(f, length, n, GLdegree):
    from basis import Poisson
    P=Poisson(n,length)
    with pytest.raises(AttributeError):
        P.u(P.nodes[0][0],P.nodes[0][1])        
    
