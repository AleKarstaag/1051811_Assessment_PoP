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
   

    

])

def test_L2norm_Helmholtz(cells_value, length, origin, f, u):
    from basis import Helmholtz 
    H=Helmholtz(cells_value, length, origin, f, u)
    cells_value2=(cells_value-2)*2+3 
    cells_value3=(cells_value2-2)*2+3
    H2=Helmholtz(cells_value2,length,origin,f,u)
    H3=Helmholtz(cells_value3,length,origin,f,u)
    e1=H.error()
    e2=H2.error()
    e3=H3.error()
    npt.assert_equal( (1/2)*H.h, H2.h )
    npt.assert_equal( (1/2)*H2.h,H3.h)
    assert e2 <= (1/4)*e1 or abs(e2-(1/4)*e1)<=1e-1
    assert e3 <= (1/4)*e2 or abs(e3-(1/4)*e2)<=2e-1
    
