"""Quadrature rules for approximating integral over
   square or triangular domains."""
import numpy as np
from scipy.special.orthogonal import p_roots


def triangle_quadrature_rule(dataframe, function,
                             vertex_1, vertex_2, vertex_3,
                             base_lenght=1):
    """Approximate double-integral of function over a
       rectangular-equilateral-triangular. """
    vertices = np.transpose(np.array([vertex_1, vertex_2, vertex_3]))
    integral = 0
    for i in range(len(dataframe)):
        barycentric_coordinates = np.array([
            float(dataframe['X'][i]),
            float(dataframe['Y'][i]),
            float(dataframe['Z'][i])])
        integral += base_lenght**2/2*function(vertices.dot(
            barycentric_coordinates))*float(dataframe['Weight'][i])
    return integral


def GaussLegendre(f, a, b, c, d, n):
    """Approximate the double-integral of the function f over the square or rectangular
       described by a<x<b and c<y<d. The degree of the Legendre polynomial is given by n.
       """
    [xi, w] = p_roots(n+1)
    x = 1/2*(1-xi)*a+1/2*(1+xi)*b
    y = 1/2*(1-xi)*c+1/2*(1+xi)*d
    G = np.zeros(((n+1), (n+1)))
    for i in range(n+1):
        for j in range(n+1):
            G[i, j] = w[i]*w[j]*f(x[i], y[j])
    G = np.sum(G)*(b-a)/2*(d-c)/2
    return G


def L2error(uh, u, a, b, c, d, n):
    def error(x, y):
        return (uh(x, y)-u(x, y))**2
    norm = np.sqrt(GaussLegendre(error, a, b, c, d, n))
    return norm

    







