import numpy as np 
import numbers
from basis.Dirichlet import nodal_basis, nodal_basis_x, nodal_basis_y
from basis.quadrature import GaussLegendre1, triangle_quadrature_rule
from numpy.core.records import array
from numpy.lib.twodim_base import tri
from scipy import linalg
from tqdm import trange
import pandas as pd


class Elliptic:


    def __init__(self, nodal_value=12, length=100,dataframe_file_name=16):
        # Number of desiderd nodes = (nodal_value-2)^2
        if not isinstance(nodal_value, numbers.Integral) and nodal_value<=2:
            raise TypeError(f"{n} is not an integer greater then 2")
        xd, yd = np.meshgrid(
            np.linspace(0, length, nodal_value), np.linspace(0, length, nodal_value))
        self.nodes = np.transpose(
            (np.concatenate(xd[1:-1, 1:-1]), np.concatenate(yd[1:-1, 1:-1])))
        self.h = xd[0, 1]
        self.a, self.b , self.c, self.d =0,length,0,length
        self.phi = lambda x, y, k: nodal_basis(x, y, self.nodes[k], self.h)
        self.phi_x = lambda x, y, k: nodal_basis_x(x, y, self.nodes[k], self.h)
        self.phi_y = lambda x, y, k: nodal_basis_y(x, y, self.nodes[k], self.h)
        self.A = "Need to run the _A() method to get numerical value."
        self.L = "Need to run the _L() method to get numerical value."
        self.dataframe=pd.read_csv(f"data/{dataframe_file_name}.csv")

    def isonthe(self,i,j):
        if self.nodes[i][0]==self.nodes[j][0]+self.h and self.nodes[i][1]==self.nodes[j][1]:
            return 'Right'
        elif self.nodes[i][0]==self.nodes[j][0]-self.h and self.nodes[i][1]==self.nodes[j][1]:
            return 'Left'
        elif self.nodes[i][0]==self.nodes[j][0] and self.nodes[i][1]==self.nodes[j][1]+self.h:
            return 'Top'
        elif self.nodes[i][0]==self.nodes[j][0] and self.nodes[i][1]==self.nodes[j][1]-self.h:
            return 'Low'
        elif self.nodes[i][0]==self.nodes[j][0] and self.nodes[i][1]==self.nodes[j][1]:
            return 'Same'
        else:
            return 'None'

    # def support(self,i,j):
    #     a=self.nodes[i][0]
    #     b=self.nodes[j][0]
    #     c=self.nodes[i][1]
    #     d=self.nodes[j][1]

    #     if (a+self.h==b and c==d) or (a - self.h == + b and c==d ):
    #         return min(a,b)-self.h/2,max(a,b)+self.h/2,c-self.h,c+self.h
    #     elif (c+self.h==d and a==b) or (c-self.h==d and a==b):
    #         return a-self.h,b+self.h,min(c,d)-self.h/2,max(c,d)+self.h/2
    #     elif a==b and c==d:
    #         return a-self.h,a+self.h,c-self.h,c+self.h
    #     else:
    #         return 0,0,0,0
class Poisson(Elliptic):


    def __init__(self, nodal_value=12, length=100):
        super().__init__(nodal_value, length)
        self.integrand_bilinear_form = lambda x, y, i, j: (
        self.phi_x(x,y,i) * self.phi_x(x,y,j) 
        + self.phi_y(x,y,i) * self.phi_y(x,y,j)
        )
    
    def _l(self, k, f = lambda x,y: np.sin(x+y), GL_degree=100):
        """Linear form of the variational formulation of 2D-Poisson PDE."""
        integrand_linear_form=lambda x,y: self.phi(x,y,k)*f(x,y)
        return GaussLegendre1(integrand_linear_form,self.a,self.b,self.c,self.d,GL_degree)

    # def _lk(self, k, f = lambda x,y: np.sin(x+y), GL_degree=100):
    #     """Linear form of the variational formulation of 2D-Poisson PDE."""
    #     integrand_linear_form=lambda x: self.phi(x[0],x[1],k)*f(x[0],x[1])
    #     return GaussLegendre1(integrand_linear_form,self.a,self.b,self.c,self.d,GL_degree)

    def _L(self, f = lambda x,y: np.sin(x+y), GL_degree = 100):
        """Returns and store as an attribute the b vector of the final linear system."""
        L=np.zeros(len(self.nodes))
        for i in range(len(self.nodes)):
            L[i]=self._l(i,f,GL_degree)
        self.L=L
        return L

    def _a(self,i,j,GL_degree=100):
        """Bilinear form of the variational formulation of 
        2D-Poisson PDE with Dirichelt BCs."""
        integrand=lambda x,y: self.integrand_bilinear_form(x,y,i,j)
        return GaussLegendre1(
            integrand,self.support(i,j),GL_degree)

    def _aT(self,i,j):
        integrand=lambda x: self.integrand_bilinear_form(x[0],x[1],i,j)
        h_shift=np.array([self.h,0])
        v_shit=np.array([0,self.h])
        if self.isonthe(i,j)=='Right':
            I=0
            I+=triangle_quadrature_rule(self.dataframe,integrand,
                self.nodes[i],self.nodes[j],self.nodes[j]+v_shit,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,
                self.nodes[i],self.nodes[j],self.nodes[i]-v_shit,self.h)
            return I
        elif self.isonthe(i,j)=='Left':
            I=0 
            I+=triangle_quadrature_rule(self.dataframe,integrand,
                self.nodes[i],self.nodes[j],self.nodes[i]+v_shit,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,
                self.nodes[i],self.nodes[j],self.nodes[j]-v_shit,self.h)
            return I
        elif self.isonthe(i,j)== 'Top':
            I=0 
            I+=triangle_quadrature_rule(self.dataframe,integrand,
                self.nodes[i],self.nodes[j],self.nodes[i]-h_shift,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,
                self.nodes[i],self.nodes[j],self.nodes[j]+h_shift,self.h)
            return I
        elif self.isonthe(i,j)=='Low':
            I=0 
            I+=triangle_quadrature_rule(self.dataframe,integrand,
                self.nodes[i],self.nodes[j],self.nodes[j]-h_shift,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,
                self.nodes[i],self.nodes[j],self.nodes[i]+h_shift,self.h)
            return I
        elif self.isonthe(i,j)=='Same':
            I=0
            I+=triangle_quadrature_rule(self.dataframe,integrand,# 1
                self.nodes[i],self.nodes[i]+h_shift,self.nodes[i]+v_shit,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,# 2
                self.nodes[i],self.nodes[i]+v_shit,self.nodes[i]+v_shit-h_shift,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,# 3
                self.nodes[i],self.nodes[i]+v_shit-h_shift,self.nodes[i]-h_shift,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,# 4
                self.nodes[i],self.nodes[i]-h_shift,self.nodes[i]-v_shit,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,# 5
                self.nodes[i],self.nodes[i]-v_shit,self.nodes[i]-v_shit+h_shift,self.h)
            I+=triangle_quadrature_rule(self.dataframe,integrand,# 6
                self.nodes[i],self.nodes[i]-v_shit+h_shift,self.nodes[i]+h_shift,self.h)
            return I
        else:
            return 0


    def _A(self,GL_degree=30):
        """Returns and store as an attribute the A matrix of the final linear system."""
        self.A=np.zeros((len(self.nodes),len(self.nodes)))
        for i in trange(len(self.nodes)): # the code uses trange to get loading bar
            for j in range(len(self.nodes)):
                self.A[i,j]=self._a(i,j,GL_degree)
        return self.A

    def _AT(self):
        """Returns and store as an attribute the A matrix of the final linear system."""
        self.AT=np.zeros((len(self.nodes),len(self.nodes)))
        for i in trange(len(self.nodes)): # the code uses trange to get loading bar
            for j in range(len(self.nodes)):
                self.AT[i,j]=self._aT(i,j)
        return self.AT

    def _U(self): # requires running self._A() and self._L() beforehand
        """Returns and store as an attribute the U vector of coefficients."""
        if isinstance(self.A,str) or isinstance(self.L,str):
            raise AttributeError(f"{type(self).__name__} object has no numerical attribute A or L."
            " Make sure you run _L() and _A() before running _U() in order to store the relevant attributes.")
        self.U=linalg.solve(self.A,self.L)
        return self.U
    
    def u(self,x,y): # requires running self._U() beforehand
        """Approximated solution u(x,y)."""    
        if isinstance(self.U,str):
            raise AttributeError(f"{type(self).__name__} object has no attribute U."
            " Make sure you run _U() before running u() in order to store the relevant attribute.")
        res=np.zeros(len(self.nodes))
        for i in trange(len(self.nodes)):
            res[i]=nodal_basis(x,y,self.nodes[i],self.h)
        return np.dot(self.U,res)
    

class Helmotz(Poisson):
    

    def __init__(self, nodal_value=12, length=100, c=lambda x,y: 4):
        super().__init__(nodal_value, length)
        self.ksq = c
        self.integrand_bilinear_form = lambda x,y,i,j: (
        self.phi_x(x, y, i) * self.phi_x(x, y, j) 
        + self.phi_y(x, y, i) * self.phi_y(x, y, j) 
        + self.phi(x, y, i) * self.ksq(x, y) )



