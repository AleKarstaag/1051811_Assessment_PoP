"""Helmotz 2D with Dirichlet BCs on square domain.

The script contains a class for solving the 2D Helmotz equation
with Dirichlet BCs on a square domain with arbitrary
location on the cartesian plane.
Examples on how to properly run the classes on
a python interpreter could be found on the report or on the test files.
"""
import numpy as np
from basis.quadrature import triangle_quadrature_rule as quad
from basis.poisson import Poisson

        
class Helmholtz (Poisson):
    """Helmholt subclass of Poisson for solving Helmholtz equation.

        Contains methods for solving Helmholtz equation on 2D Square Domain
        by finite element method with piecewise linear functions.

        Parameteres
        ----------
        nodal_value : integer
            The desired number of nodes in the first row of nodes;
        length : float
            The desired lenght of the side of square domain.
        origin : 2-sized list
            The desired point from where the square-domain will start.
        f: function(x : float, y: float) -> float
            The f(x,y) function of the Poisson equation.
        u: function(x : float, y: float) -> float (Optional)
            The u(x,y) analytical solution of the Poisson equation
        c: function(x : float, y: float) -> float

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
            Barycentric coordinates and weights for
            16 points triangular quadrature.
        data46 : pd.core.frame.DataFrame
            Barycentric coordinates and weights for
            46 points trinagular quadrature.
        domain : list( float, float, float, float )
            Contains information descibing the location of the square-domain.
        A : str or np.ndarray
            Stiffness matrix
        b : str or np.ndarray
            Right hand side of linear system with stiffness matrix.
        U : str or np.ndarray
            Solution array of linear system.
        L2error : float (Optional)
            L2norm of the error between analytical
            solution and Galerkin approximation.
        integrand_bilinear_form: function(x : float, y : float ,
                                          i : int, j : int) -> float
            The integrand of the bilenar form at nodes i and j.
        f = function (x : float, y : float) -> float
            The f(x,y) function of the Poisson equation.
        u = function (x : float, y : float) -> float (Optional)
            The u(x,y) analytical solution of the Poisson equation.
        c = function (x : float, y : float) -> float

        """

    def __init__(self, nodal_value=12, length=100, origin=[0, 0],
                 f=lambda x, y:
                 (2*np.pi**2/(100**2)+1) * np.sin(x*np.pi/100)
                 * np.sin(y*np.pi/100),
                 u=lambda x, y: np.sin(x*np.pi/100) * np.sin(y*np.pi/100),
                 c=1):
        super().__init__(nodal_value, length, origin, f, u)
        self.c = c
        self.integrand_bilinear_form = lambda x, y, i, j: (
                self.phi_x(x, y, i) * self.phi_x(x, y, j)
                + self.phi_y(x, y, i) * self.phi_y(x, y, j)
                + self.phi(x, y, i) * self.phi(x, y, j) * self.c)

    def _a(self, i, j):
        """Approximate the bilinear form only on the correspoding support."""
        def integrand(x):
            return self.integrand_bilinear_form(x[0], x[1], i, j)
        h_shift = np.array([self.h, 0])
        v_shift = np.array([0, self.h])
        if self.isonthe(i, j) == 'Right':
            integral = 0
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[j]+v_shift, self.h)
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[i]-v_shift, self.h)
            return integral
        elif self.isonthe(i, j) == 'Left':
            integral = 0
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[i]+v_shift, self.h)
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[j]-v_shift, self.h)
            return integral
        elif self.isonthe(i, j) == 'Top':
            integral = 0
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[i]-h_shift, self.h)
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[j]+h_shift, self.h)
            return integral
        elif self.isonthe(i, j) == 'Low':
            integral = 0
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[j]-h_shift, self.h)
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[i]+h_shift, self.h)
            return integral
        elif self.isonthe(i, j) == 'TopLeft':
            integral = 0
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[i]+h_shift, self.h)
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[j]-h_shift, self.h)
            return integral
        elif self.isonthe(i, j) == 'LowRight':
            integral = 0
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[i]+v_shift, self.h)
            integral += quad(self.data46, integrand,
                             self.nodes[i], self.nodes[j],
                             self.nodes[j]-v_shift, self.h)
            return integral
        elif self.isonthe(i, j) == 'Same':
            integral = 0
            integral += quad(self.data46, integrand,  # 1
                             self.nodes[i], self.nodes[i]+h_shift,
                             self.nodes[i]+v_shift, self.h)
            integral += quad(self.data46, integrand,  # 2
                             self.nodes[i], self.nodes[i]+v_shift,
                             self.nodes[i]+v_shift-h_shift, self.h)
            integral += quad(self.data46, integrand,  # 3
                             self.nodes[i], self.nodes[i]+v_shift-h_shift,
                             self.nodes[i]-h_shift, self.h)
            integral += quad(self.data46, integrand,  # 4
                             self.nodes[i], self.nodes[i]-h_shift,
                             self.nodes[i]-v_shift, self.h)
            integral += quad(self.data46, integrand,  # 5
                             self.nodes[i], self.nodes[i]-v_shift,
                             self.nodes[i]-v_shift+h_shift, self.h)
            integral += quad(self.data46, integrand,  # 6
                             self.nodes[i], self.nodes[i]-v_shift+h_shift,
                             self.nodes[i]+h_shift, self.h)
            return integral
        else:
            return 0
