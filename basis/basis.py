import numpy as np 
import numbers
from scipy.integrate import quad
from scipy.special.orthogonal import p_roots

def dim(matrix):
    """Returns a tuple identifying the number of rows and number of columns of the array."""
    if not type(matrix)==np.ndarray:
        raise NotImplementedError('Make sure the argument is a numpy.ndarray')
    else:
        return (len(matrix),len(matrix[0]))
    
def nodal_basis(x,y,nodal_point,h):
    """Finite basis function."""
    x_i=nodal_point[0]
    y_j=nodal_point[1]
    if (y-y_j)>=0 and (x-x_i)>=0 and (y-y_j) <= -(x-x_i)+h: #1
        return 1 - (x-x_i)/h - (y-y_j)/h
    elif -h <= (x-x_i) <= 0 and (y-y_j) >= -(x-x_i) and h>= (y-y_j) >= 0: #2
        return 1 - (y-y_j)/h
    elif -h <= (x-x_i) <= 0 and (y-y_j) <= -(x-x_i) and h>=(y-y_j)>=0: #3
        return 1 - (x_i - x)/h
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h : #4  
        return 1 - (x_i-x)/h - (y_j-y)/h
    elif (h+x_i)>=x>= x_i and (y-y_j) <= -(x-x_i) and (-h+y_j)<=y<=y_j: #5
        return 1 - (y_j - y)/h
    elif (h+x_i)>=x>=x_i and (y-y_j) >= -(x-x_i) and (-h+y_j) <=y<=y_j: #6
        return 1 - (x-x_i)/h
    elif x==x_i and y==y_j:
        return 1    
    else:
        return 0

def nodal_basis_codomain(mesh,nodal_point,h):
    xd=mesh[0]
    yd=mesh[1]
    res=np.zeros(dim(xd))
    for i in range(dim(xd)[0]):
        for j in range(dim(xd)[1]):  
            res[i,j]= nodal_basis(xd[i,j],yd[i,j],nodal_point,h)
    return res

def nodal_basis_y(x,y,nodal_point,h):
    """Finite basis function."""
    x_i=nodal_point[0]
    y_j=nodal_point[1]
    if (y-y_j)>=0 and (x-x_i)>=0 and (y-y_j) <= -(x-x_i)+h: #1
        return -1/h
    elif -h <= (x-x_i) <= 0 and (y-y_j) >= -(x-x_i) and h>= (y-y_j) >= 0: #2
        return -1/h 
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h : #4  
        return 1/h
    elif (h+x_i)>=x>= x_i and (y-y_j) <= -(x-x_i) and (-h+y_j)<=y<=y_j: #5
        return 1/h
    else:
        return 0

def nodal_basis_x(x,y,nodal_point,h):
    """Finite basis function."""
    x_i=nodal_point[0]
    y_j=nodal_point[1]
    if (y-y_j)>=0 and (x-x_i)>=0 and (y-y_j) <= -(x-x_i)+h: #1
        return -1/h
    elif -h <= (x-x_i) <= 0 and (y-y_j) <= -(x-x_i) and h>=(y-y_j)>=0: #3
        return 1/h
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h : #4  
        return 1/h
    elif (h+x_i)>=x>=x_i and (y-y_j) >= -(x-x_i) and (-h+y_j) <=y<=y_j: #6
        return -1/h
    else:
        return 0

def GaussLegendre(phi,f,a,b,c,d,n=7):
    [xi,w]=p_roots(n+1)
    x=1/2*(1-xi)*a+1/2*(1+xi)*b
    y=1/2*(1-xi)*c+1/2*(1+xi)*d
    G=np.zeros(((n+1),(n+1)))
    for i in range(n+1):
        for j in range(n+1):
            G[i,j]= w[i]*w[j]*phi(x[i],y[j])*f(x[i],y[j])
    G=np.sum(G)*(b-a)/2*(d-c)/2
    return G
