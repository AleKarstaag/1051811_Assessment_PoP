import numpy as np 
import numbers
from basis.basis import basis, GaussLegendre, finite_basis
from scipy.special.orthogonal import p_roots

class SqGrid:


    def __init__(self,l=100,n=4,m=25):
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

   
    def _evaluate(self,k):
        """"Returns a grid representing the codomain of the piecewise-continuos basis
        linear function corresponding to the k-th nodal point.
        """
        res=basis(self.mesh,self.nodes[k],self.h)
        return res

    
    def _check(self,k):
        """Continuos piecewise linear basis functions should map their corresponding nodes to
        the value: 1. The computational grid doesn't necesessarily generate the nodal point.
        This function has been created only to check approximately whether the function has been coded correctly.
        """ 
        res = self.evaluate(k)
        if 0.965<=res[round(self.nodes[k][1]),round(self.nodes[k][0])]<=1:
            return f"({round(self.nodes[k][0])},{round(self.nodes[k][1])}) seems correct."
        elif 0.94<=res[round(self.nodes[k][1]),round(self.nodes[k][0])]<0.965:
            raise Warning(f"({round(self.nodes[k][0])},{round(self.nodes[k][1])}) might not be correct.")
        else:
            raise EstimationError("Not correct.")

    def _machine_errors(self,k):
        """Returns coordinate of points where the function 
           finite_basis maps to a negative value.
           This has been used often to 
           check whether the function has been defined correctly.
           """
        res=self.evaluate(k)
        logic=[]
        for i in range(len(res)):
            for j in range(len(res[0])):
                if res[i,j]<0:
                    logic.append(f"({i},{j})")
        return np.transpose(logic)

    # def _clean_evaluate(self,k):
    #     res=self.evaluate(k)
    #     for i in range(len(res)):
    #         for j in range(len(res[0])):
    #             if res[i,j]<0:
    #                 res[i,j]=0
    #     return res

    """Integrators"""

    def _inner(self,k,n=7,f=lambda x,y: 1):
        """Gauss-Legendre approximation (n degree) of the integral of
           the k-th piecewise continuos linear function times 
           an arbitrary function f over the whole domain of the problem.
           In short: 
           let ( u , v )= int_{Omega} u v
           The method approximates:
                    ( phi_k(x,y) , f(x,y) ) if 3 inputs are given
                    ||phi_k(x,y)||^2 if 1 or 2 inputs are given
                    """
        a,b=self.mesh[0][0,0],self.mesh[0][0,-1]
        c,d=self.mesh[1][0,0],self.mesh[1][-1,0]
        phi=lambda x,y: finite_basis(x,y,self.nodes[k],self.h)
        G=GaussLegendre(phi,f,a,b,c,d,n)
        return G

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
    
class EstimationError(Exception):
    pass


    

    




  
    

