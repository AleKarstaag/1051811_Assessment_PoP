import numpy as np
from scipy.special.orthogonal import p_roots

def GaussLegendre1(function,a,b,c,d,n=7):
    [xi,w]=p_roots(n+1)
    x=1/2*(1-xi)*a+1/2*(1+xi)*b
    y=1/2*(1-xi)*c+1/2*(1+xi)*d
    G=np.zeros(((n+1),(n+1)))
    for i in range(n+1):
        for j in range(n+1):
            G[i,j]= w[i]*w[j]*function(x[i],y[j])
    return np.sum(G)*(b-a)/2*(d-c)/2
