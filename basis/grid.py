import numpy as np 
import numbers
from basis.Dirichlet import nodal_basis, nodal_basis_x, nodal_basis_y
from basis.quadrature import GaussLegendre
from scipy import linalg
class Elliptic:
    def __init__(self,n=4,l=100,m=25,f=None):
        if not isinstance(n,numbers.Integral):
            raise TypeError(f"{n} is not an integer")
        xd , yd = np.meshgrid(np.linspace(0,l,n),np.linspace(0,l,n))
        nodes = np.transpose((np.concatenate(xd[1:-1,1:-1]),np.concatenate(yd[1:-1,1:-1])))
        #notice: Number of nodes = (n-2)(n-2)
        xv,yv = np.meshgrid(np.linspace(0,l,round(n*m)),np.linspace(0,l,round(n*m)))
        self.h = xd[0,1]
        self.dim = (len(xv),len(yv[0]))
        self.nodes = nodes
        self.mesh = xv,yv
        self.area= xv[0,-1]*yv[-1,0]
        self.a,self.b=self.mesh[0][0,0],self.mesh[0][0,-1]
        self.c,self.d=self.mesh[1][0,0],self.mesh[1][-1,0]
        self.f=f

class Poisson(Elliptic):

   
    def _F(self,k,f=None,n=20):
        """Linear form of the variational formulation of 2D-Poisson PDE."""
        if f==None:
            f=self.f     
        phi=lambda x,y: nodal_basis(x,y,self.nodes[k],self.h)
        G=GaussLegendre(phi,f,self.a,self.b,self.c,self.d,n)
        return G

    def _b(self,f=lambda x,y: 1,n=100):
        """Returns the b vector of the final linear system."""
        b=np.zeros(len(self.nodes))
        for i in range(len(self.nodes)):
            b[i]=self._F(i,f,n)
        return b

    def _a(self,i,j,n=100):
        """Bilinear form of the variational formulation of 
        2D-Poisson PDE with Dirichelt BCs."""

        phi_i_x=lambda x,y: nodal_basis_x(x,y,self.nodes[i],self.h)
        phi_j_x=lambda x,y: nodal_basis_x(x,y,self.nodes[j],self.h)
        phi_i_y=lambda x,y: nodal_basis_y(x,y,self.nodes[i],self.h)
        phi_j_y=lambda x,y: nodal_basis_y(x,y,self.nodes[j],self.h)
        Gx=GaussLegendre(phi_i_x,phi_j_x,self.a,self.b,self.c,self.d,n)
        Gy=GaussLegendre(phi_i_y,phi_j_y,self.a,self.b,self.c,self.d,n)
        return Gx+Gy
    
    def _A(self,n=300):
        """Returns the A matrix of the final linear system."""
        A=np.zeros((len(self.nodes),len(self.nodes)))
        for i in range(len(self.nodes)):
            for j in range(len(self.nodes)):
                A[i,j]=self._a(i,j,n)
        return A

    def _U(self,f=lambda x,y: 1,n=300):
        U=linalg.solve(self._A(n),self._b(f,round(n/3)))
        return U
    
    def _u(self,x,y,f,n):
        U=self._U(f,n)
        res=np.zeros(len(self.nodes))
        for i in range(len(self.nodes)):
            res[i]=nodal_basis(x,y,self.nodes[i],self.h)
        return np.dot(U,res)

class Helmotz(Elliptic):
    pass
