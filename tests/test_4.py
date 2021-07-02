from unicodedata import decimal
import pytest
import numpy as np
import numpy.testing as npt
"""Test Galerkin approximation of the Helmotz equation."""

@pytest.mark.parametrize("cells_value, length, origin,f,u",[
    (5,100,[0,0],
    lambda x,y: (2*np.pi**2/(100**2)+1) * np.sin(x*np.pi/100)*np.sin(y*np.pi/100),
    lambda x,y: np.sin(x*np.pi/100) * np.sin(y*np.pi/100)),

    (5,4,[-2,-2],
    lambda x,y: -2*(y**2-4)-2*(x**2-4)+(x**2-4)*(y**2-4),
    lambda x,y: (x**2-4)*(y**2-4)),

    (5,2,[-1,-1],
    lambda x,y: -(2*np.exp(x**2)+4*np.exp(x**2)*x**2)*(y**2-1)
                -2*(np.exp(x**2)-np.e)
                +(np.exp(x**2)-np.e)*(y**2-1),
    lambda x,y: (np.exp(x**2)-np.e)*(y**2-1))
    # (5,4,[-2,-2],
    # lambda x,y: -2*(y**2-4)-2*(x**2-4),
    # lambda x,y: (x**2-4)*(y**2-4)),

    # (5,2,[-1,-1],
    # lambda x,y: -(2*np.exp(x**2)+4*np.exp(x**2)*x**2)*(y**2-1)-2*(np.exp(x**2)-np.e),
    # lambda x,y: (np.exp(x**2)-np.e)*(y**2-1))

    

])

def test_L2norm(cells_value, length, origin, f, u):
    from basis import Helmotz
    H=Helmotz(cells_value, length, origin, f, u)
    cells_value_new=(cells_value-2)*2+3 
    cells_value_new_new=(cells_value_new-2)*2+3
    H_new=Helmotz(cells_value_new,length,origin,f,u)
    H_new_new=Helmotz(cells_value_new_new,length,origin,f,u)
    error=H.error()
    error_new=H_new.error()
    error_new_new=H_new_new.error()
    npt.assert_equal( (1/2)*H.h, H_new.h )
    npt.assert_equal( (1/2)*H_new.h,H_new_new.h)
    npt.assert_almost_equal((1/4)*error, error_new, decimal=0)
    npt.assert_almost_equal((1/4)*error_new, error_new_new, decimal=0)
    
