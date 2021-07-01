from unicodedata import decimal
import pytest
import numpy as np
import numpy.testing as npt
"""Test Galerkin approximation of solution of the Poisson equation."""

@pytest.mark.parametrize("cells_value, length, origin,f,u",[
    (5,100,[0,0],
    lambda x,y: (2*np.pi**2/(100**2)) * np.sin(x*np.pi/100)*np.sin(y*np.pi/100),
    lambda x,y: np.sin(x*np.pi/100) * np.sin(y*np.pi/100)),

    (9,4,[-2,-2],lambda x,y: -2*(y**2-4)-2*(x**2-4),lambda x,y: (x**2-4)*(y**2-4)),

    (9,2,[-1,-1],
    lambda x,y: -(2*np.exp(x**2)+4*np.exp(x**2)*x**2)*(y**2-1)-2*(np.exp(x**2)-np.e),
    lambda x,y: (np.exp(x**2)-np.e)*(y**2-1))

])

def test_L2norm(cells_value, length, origin, f, u):
    from basis import Poisson
    P=Poisson(cells_value, length, origin, f, u)
    cells_value_new=(cells_value-2)*2+3 
    P_new=Poisson(cells_value_new,length,origin,f,u)
    error=P.error()
    error_new=P_new.error()
    npt.assert_equal( (1/2)*P.h, P_new.h )
    npt.assert_almost_equal((1/4)*error, error_new, decimal=1)
    
