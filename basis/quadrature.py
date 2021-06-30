import numpy as np
from scipy.special.orthogonal import p_roots
import pandas as pd

def GaussLegendre1(function,interval,n=7):
    [xi,w]=p_roots(n+1)
    x=1/2*(1-xi)*interval[0]+1/2*(1+xi)*interval[1]
    y=1/2*(1-xi)*interval[2]+1/2*(1+xi)*interval[3]
    G=np.zeros(((n+1),(n+1)))
    for i in range(n+1):
        for j in range(n+1):
            G[i,j]= w[i]*w[j]*function(x[i],y[j])
    return np.sum(G)*(interval[1]-interval[0])/2*(interval[3]-interval[2])/2

def triangle_quadrature_rule(dataframe,function,vertex_1,vertex_2,vertex_3,base_lenght=1,name_file=16):
    vertices=np.transpose(np.array([vertex_1,vertex_2,vertex_3]))
    # quadrature_nodes=[]
    # weights=[]
    I=0
    for i in range(len(dataframe)):
        barycentric_coordinates=np.array([
            float(dataframe['X'][i]),
            float(dataframe['Y'][i]),
            float(dataframe['Z'][i])])
        # quadrature_nodes.append(vertices.dot(barycentric_coordinates))
        # weights.append(float(dataframe['Weight'][i]))
        I += base_lenght**2/2*function(vertices.dot(barycentric_coordinates))*float(dataframe['Weight'][i])
    
    return I

    







