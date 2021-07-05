"""Elliptic 2D with Dirichlet BCs on square domain.

The script contains a class for initialsing an
elliptic PDE solver on a 2D square domain with
Dirichlet boundaries conditions .
"""
import numpy as np
from basis.Dirichlet import nodal_basis, nodal_basis_x, nodal_basis_y
import pandas as pd


class Elliptic:
    """Initialise data for solving an Elliptic PDE in a 2D square domain.

    Contains information about the domain and the cells that will
    be used in the finite element approximation of the domain.

    Parameteres
    ----------
    nodal_value : integer
        The desired number of nodes in the first row of nodes;
    length : float
        The desired lenght of the side of square domain.
    origin : 2-sized list
        The desired point from where the square-domain will start,


    Attributes
    ----------
    nodes : np.ndarray
        Cartesian coordinates of the interior nodes.
    h : float
        Diameter
    length : float
        Length of the side of the square-domain
    phi : function(x : float , y : float, k : integer) -> float
        Piecewise linear function of the finite element method
    phi_x : function(x : float , y : float, k : integer) -> float
        Partial derivative w.r.t x of phi
    phi_y : function(x : float , y : float, k : integer) -> float
        Partial derivative w.r.t y of phi
    data16 : pd.core.frame.DataFrame
        Barycentric coordinates and weights for 16 points triangular quadrature
    data46 : pd.core.frame.DataFrame
        Barycentric coordinates and weights for 46 points trinagular quadrature
    domain : list( float, float, float, float )
        Contains information descibing the location of the square-domain.

    """

    def __init__(self, nodal_value=12, length=100, origin=[0, 0]):
        xd, yd = np.meshgrid(
            np.linspace(origin[0], length+origin[0], nodal_value),
            np.linspace(origin[1], length+origin[1], nodal_value))
        self.nodes = np.transpose(
            (np.concatenate(xd[1:-1, 1:-1]), np.concatenate(yd[1:-1, 1:-1])))
        self.h = xd[0, 1]-xd[0, 0]
        self.length = length
        self.phi = lambda x, y, k: nodal_basis(x, y, self.nodes[k], self.h)
        self.phi_x = lambda x, y, k: nodal_basis_x(x, y, self.nodes[k], self.h)
        self.phi_y = lambda x, y, k: nodal_basis_y(x, y, self.nodes[k], self.h)
        self.data16 = pd.read_csv("data/16.csv")
        self.data46 = pd.read_csv("data/46.csv")
        self.domain = [origin[0], length+origin[0],
                       origin[1], length+origin[1]]

    def __repr__(self):
        return self.__class__.__name__ + " with x in " + repr(
            self.domain[0]) + " , " + " y in "+repr(
            self.domain[1]) + " and " + repr(len(self.nodes)) + " cells "

    def isonthe(self, i, j, tol=1e-10):
        """Returns a string saying where the i-th interior node is with
        respect to the j-th one."""

        if not -len(self.nodes) < i < len(self.nodes):
            raise ValueError(f"{i} is not between {-len(self.nodes)}"
                             f"and {len(self.nodes)}")
        if not -len(self.nodes) < j < len(self.nodes):
            raise ValueError(f"{j} is not between {-len(self.nodes)}"
                             f" and {len(self.nodes)}")
        if self.nodes[i][1] == self.nodes[j][1]:
            if abs(self.nodes[i][0]-self.nodes[j][0] - self.h) <= tol:
                return 'Right'
            elif abs(self.nodes[i][0]-self.nodes[j][0]+self.h) <= tol:
                return 'Left'
            elif self.nodes[i][0] == self.nodes[j][0]:
                return 'Same'
            else:
                return 'None'
        elif self.nodes[i][0] == self.nodes[j][0]:
            if abs(self.nodes[i][1]-self.nodes[j][1] - self.h) <= tol:
                return 'Top'
            elif abs(self.nodes[i][1]-self.nodes[j][1]+self.h) <= tol:
                return 'Low'
            else:
                return 'None'
        elif abs(self.nodes[i][1] - self.nodes[j][1] - self.h) <= tol:
            if abs(self.nodes[i][0]-self.nodes[j][0] + self.h) <= tol:
                return 'TopLeft'
            else:
                return 'None'
        elif abs(self.nodes[i][1] - self.nodes[j][1] + self.h) <= tol:
            if abs(self.nodes[i][0]-self.nodes[j][0] - self.h) <= tol:
                return 'LowRight'
            else:
                return 'None'
        else:
            return 'None'

