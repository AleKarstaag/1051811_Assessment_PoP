import pytest
import numpy as np
import numpy.testing as npt

@pytest.mark.parametrize("function, vertex_1, vertex_2, vertex_3,ans", [
    (lambda x: x[0]**2,[0,0],[1,0],[1,1],1/4),
    (lambda x: x[0]**2,[0,0],[1,0],[0,1],-1/4+1/3),
    (lambda x: np.sin(x[0]+x[1]),[0,0],[1,0],[0,1],-np.cos(1)+np.sin(1)),
    (lambda x: x[0]**2,[-1,0],[0,-1],[-1,-1],1/4),
    (lambda x: 1,[0,0],[1,0],[0,1],1/2),
    (lambda x: np.exp(x[0]),[0,0],[1,0],[1,1],1),
    (lambda x: 2*np.exp(x[0])*x[1],[0,0],[1,0],[1,1],np.exp(1)-2)
    ])
    
def test_isontheright(function, vertex_1,vertex_2,vertex_3,ans):
    """Test """
    from basis import EllipticDirichlet
    I,X,W=triangle_quadrature_rule(function,vertex_1,vertex_2,vertex_3)
    npt.assert_equal(I,ans , decimal=7)