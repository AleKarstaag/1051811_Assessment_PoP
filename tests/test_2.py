import pytest
import numpy as np
import numpy.testing as npt
import pandas as pd

@pytest.mark.parametrize("dataframe_file_name,"
    "function, vertex_1, vertex_2, vertex_3,ans", [
    (16,lambda x: x[0]**2,[0,0],[1,0],[1,1],1/4),
    (16,lambda x: x[0]**2,[0,0],[1,0],[0,1],-1/4+1/3),
    (16,lambda x: np.sin(x[0]+x[1]),[0,0],[1,0],[0,1],-np.cos(1)+np.sin(1)),
    (16,lambda x: x[0]**2,[-1,0],[0,-1],[-1,-1],1/4),
    (16,lambda x: 1,[0,0],[1,0],[0,1],1/2),
    (16,lambda x: np.exp(x[0]),[0,0],[1,0],[1,1],1),
    (16,lambda x: 2*np.exp(x[0])*x[1],[0,0],[1,0],[1,1],np.exp(1)-2),
    ])
    
def test_triangle_quadrature_rule16(dataframe_file_name,function, vertex_1,vertex_2,vertex_3,ans):
    """Test """
    dataframe=pd.read_csv(f"data/{dataframe_file_name}.csv")
    from basis import triangle_quadrature_rule
    I=triangle_quadrature_rule(dataframe,function,vertex_1,vertex_2,vertex_3)
    npt.assert_almost_equal(I,ans , decimal=10) 

@pytest.mark.parametrize("dataframe_file_name,"
    "function, vertex_1, vertex_2, vertex_3,ans", [
    (46,lambda x: x[0]**2,[0,0],[1,0],[1,1],1/4),
    (46,lambda x: x[0]**2,[0,0],[1,0],[0,1],-1/4+1/3),
    (46,lambda x: np.sin(x[0]+x[1]),[0,0],[1,0],[0,1],-np.cos(1)+np.sin(1)),
    (46,lambda x: x[0]**2,[-1,0],[0,-1],[-1,-1],1/4),
    (46,lambda x: 1,[0,0],[1,0],[0,1],1/2),
    (46,lambda x: np.exp(x[0]),[0,0],[1,0],[1,1],1),
    (46,lambda x: 2*np.exp(x[0])*x[1],[0,0],[1,0],[1,1],np.exp(1)-2)
    ])
def test_triangle_quadrature_rule46(dataframe_file_name,function, vertex_1,vertex_2,vertex_3,ans):
    """Test """
    dataframe=pd.read_csv(f"data/{dataframe_file_name}.csv")
    from basis import triangle_quadrature_rule
    I=triangle_quadrature_rule(dataframe,function,vertex_1,vertex_2,vertex_3)
    npt.assert_almost_equal(I,ans , decimal=15) 