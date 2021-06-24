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


    # def _evaluate(self,k):
    #     """"Returns a grid representing the codomain of the piecewise-continuos basis
    #     linear function corresponding to the k-th nodal point.
    #     """
    #     res=nodal_basis_codomain(self.mesh,self.nodes[k],self.h)
    #     return res

    # def _check(self,k):
    #     """Continuos piecewise linear basis functions should map their corresponding nodes to
    #     the value: 1. The computational grid doesn't necesessarily generate the nodal point.
    #     This function has been created only to check approximately whether the function has been coded correctly.
    #     """ 
    #     res = self.evaluate(k)
    #     if 0.965<=res[round(self.nodes[k][1]),round(self.nodes[k][0])]<=1:
    #         return f"({round(self.nodes[k][0])},{round(self.nodes[k][1])}) seems correct."
    #     elif 0.94<=res[round(self.nodes[k][1]),round(self.nodes[k][0])]<0.965:
    #         raise Warning(f"({round(self.nodes[k][0])},{round(self.nodes[k][1])}) might not be correct.")
    #     else:
    #         raise EstimationError("Not correct.")

    # def _machine_errors(self,k):
    #     """Returns coordinate of points where the function 
    #        finite_basis maps to a negative value.
    #        This has been used often to 
    #        check whether the function has been defined correctly.
    #        """
    #     res=self.evaluate(k)
    #     logic=[]
    #     for i in range(len(res)):
    #         for j in range(len(res[0])):
    #             if res[i,j]<0:
    #                 logic.append(f"({i},{j})")
    #     return np.transpose(logic)

    # def _clean_evaluate(self,k):
    #     res=self.evaluate(k)
    #     for i in range(len(res)):
    #         for j in range(len(res[0])):
    #             if res[i,j]<0:
    #                 res[i,j]=0
    #     return res
    # def _MC(self,k):
    #     """Integrate the continuos piecewise linear function corresponding to the 
    #        k-th node over the whole domain using a method similar to Monte-Carlo.
    #        It's probably the method that requires the most computational power and on the
    #        other side it probably doesn't give good results.
    #        """
    #     return np.sum(self.evaluate(k))/np.size(self.evaluate(k))*self.area

    # def _integr(self,k):
    #     """Integrate using a quadrature rule over the x-axis and then brings back to a 
    #        Monte Carlo 1-D integration. Thsi requires less computation and it's probably 
    #        more precise than the previous one.
    #        """
    #     x_i=self.nodes[k][0]
    #     y_j=self.nodes[k][1]
    #     I=[]
    #     for y in self.mesh[1][:,0]:
    #         I.append(quad(phi_ij,self.mesh[0][0,0],self.mesh[0][0,-1],args=(y,x_i,y_j,self.h))[0])
    #     return np.sum(I)/np.size(I)*self.mesh[1][-1,0]

