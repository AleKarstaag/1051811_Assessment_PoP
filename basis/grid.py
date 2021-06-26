import numpy as np 
import numbers
from basis.Dirichlet import nodal_basis, nodal_basis_x, nodal_basis_y
from basis.quadrature import GaussLegendre, GaussLegendre3
from scipy import linalg


class Elliptic:


    def __init__(self, n=4, length=100, m=25):
        if not isinstance(n, numbers.Integral):
            raise TypeError(f"{n} is not an integer")
        xd, yd = np.meshgrid(
            np.linspace(0, length, n), np.linspace(0, length, n))
        nodes = np.transpose(
            (np.concatenate(xd[1:-1, 1:-1]), np.concatenate(yd[1:-1, 1:-1])))
        # notice: Number of nodes = (n-2)(n-2)
        xv, yv = np.meshgrid(
            np.linspace(0, length, round(n*m)),
            np.linspace(0, length, round(n*m)))
        self.h = xd[0, 1]
        self.dim = (len(xv), len(yv[0]))
        self.nodes = nodes
        self.mesh = xv, yv
        self.area = xv[0, -1]*yv[-1, 0]
        self.a, self.b , self.c, self.d =xv[0,0],xv[0,-1],yv[0,0],yv[-1,0]
        
class Poisson(Elliptic):
    

    def _l(self, k, f=lambda x,y: 1, n=20):
        """Linear form of the variational formulation of 2D-Poisson PDE."""
        phi_k=lambda x,y: nodal_basis(x,y,self.nodes[k],self.h)
        G=GaussLegendre(phi_k,f,self.a,self.b,self.c,self.d,n)
        return G

    def _L(self,f=lambda x,y: 1,n=100):
        """Returns and store as an attribute the b vector of the final linear system."""
        L=np.zeros(len(self.nodes))
        for i in range(len(self.nodes)):
            L[i]=self._l(i,f,n)
        self.L=L
        return L

    def _a(self,i,j,n=100):
        """Bilinear form of the variational formulation of 
        2D-Poisson PDE with Zero BCs."""

        phi_i_x=lambda x,y: nodal_basis_x(x,y,self.nodes[i],self.h)
        phi_j_x=lambda x,y: nodal_basis_x(x,y,self.nodes[j],self.h)
        phi_i_y=lambda x,y: nodal_basis_y(x,y,self.nodes[i],self.h)
        phi_j_y=lambda x,y: nodal_basis_y(x,y,self.nodes[j],self.h)
        Gx=GaussLegendre(phi_i_x,phi_j_x,self.a,self.b,self.c,self.d,n)
        Gy=GaussLegendre(phi_i_y,phi_j_y,self.a,self.b,self.c,self.d,n)
        return Gx+Gy
    
    def _A(self,n=300):
        """Returns and store as an attribute the A matrix of the final linear system."""
        self.A=np.zeros((len(self.nodes),len(self.nodes)))
        for i in range(len(self.nodes)):
            for j in range(len(self.nodes)):
                self.A[i,j]=self._a(i,j,n)
        return self.A

    def _U(self): # requires running self._A() and self._L() beforehand
        """Returns and store as an attribute the U vector of coefficients."""
        if AttributeError:
            raise AttributeError(f"{type(self).__name__} object has no attribute A or L."
            "Make sure you run _L() and _A() before running ._U() in order to store the relevant attributes.")
        
        self.U=linalg.solve(self.A,self.L)
        return self.U
    
    def u(self,x,y): # requires running self._U() beforehand
        """Approximated solution u(x,y)."""    
        if AttributeError:
            raise AttributeError(f"{type(self).__name__} object has no attribute U."
            "Make sure you run _U() before running u() in order to store the relevant attribute.")
        res=np.zeros(len(self.nodes))
        for i in range(len(self.nodes)):
            res[i]=nodal_basis(x,y,self.nodes[i],self.h)
        return np.dot(self.U,res)
    

class Helmotz(Poisson):
    

    def __init__(self, n=4, length=100, m=25,ksq = lambda x,y: 4**2 ):
        super().__init__(n, length, m)
        self.ksq=ksq

    def _a(self,i,j,n=100):
        """Bilinear form of the variational formulation of 
        2D-Hoisson PDE with Zero BCs."""

        phi_i_x=lambda x,y: nodal_basis_x(x,y,self.nodes[i],self.h)
        phi_j_x=lambda x,y: nodal_basis_x(x,y,self.nodes[j],self.h)
        phi_i_y=lambda x,y: nodal_basis_y(x,y,self.nodes[i],self.h)
        phi_j_y=lambda x,y: nodal_basis_y(x,y,self.nodes[j],self.h)
        Gx=GaussLegendre(phi_i_x,phi_j_x,self.a,self.b,self.c,self.d,n)
        Gy=GaussLegendre(phi_i_y,phi_j_y,self.a,self.b,self.c,self.d,n)
        helmotz=GaussLegendre3(lambda x,y: nodal_basis(x,y,self.nodes[i],self.h),
                        lambda x,y: nodal_basis(x,y,self.nodes[j],self.h),
                        self.ksq,self.a,self.b,self.c,self.d,n)
        return Gx+Gy+helmotz

