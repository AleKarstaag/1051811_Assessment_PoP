import numpy as np
from scipy.special.orthogonal import p_roots
import pandas as pd

def triangle_quadrature_rule(dataframe,function,vertex_1,vertex_2,vertex_3,base_lenght=1,name_file=16):
    vertices=np.transpose(np.array([vertex_1,vertex_2,vertex_3]))
    I=0
    for i in range(len(dataframe)):
        barycentric_coordinates=np.array([
            float(dataframe['X'][i]),
            float(dataframe['Y'][i]),
            float(dataframe['Z'][i])])
        I += base_lenght**2/2*function(vertices.dot(barycentric_coordinates))*float(dataframe['Weight'][i])
    return I




    







