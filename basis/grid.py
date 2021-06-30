import numpy as np 
import numbers
from basis.Dirichlet import nodal_basis, nodal_basis_x, nodal_basis_y
from basis.quadrature import triangle_quadrature_rule as quad
from scipy import linalg
from tqdm import trange
import pandas as pd


class Elliptic:


    def __init__(self, nodal_value=12, length=100):
        # Number of desiderd nodes = (nodal_value-2)^2
        if not isinstance(nodal_value, numbers.Integral) and nodal_value<=2:
            raise TypeError(f"{n} is not an integer greater then 2")
        xd, yd = np.meshgrid(
            np.linspace(0, length, nodal_value), np.linspace(0, length, nodal_value))
        self.nodes = np.transpose(
            (np.concatenate(xd[1:-1, 1:-1]), np.concatenate(yd[1:-1, 1:-1])))
        self.h = xd[0, 1]
        self.length=length
        self.phi = lambda x, y, k: nodal_basis(x, y, self.nodes[k], self.h)
        self.phi_x = lambda x, y, k: nodal_basis_x(x, y, self.nodes[k], self.h)
        self.phi_y = lambda x, y, k: nodal_basis_y(x, y, self.nodes[k], self.h)
        self.A = "Need to run the _A() method to get numerical value."
        self.L = "Need to run the _L() method to get numerical value."
        self.U = "Need to run the _U() method to get numerical value."
        self.data16=pd.read_csv(f"data/{16}.csv")
        self.data46=pd.read_csv(f"data/{46}.csv")

    def isonthe(self,i,j,tol=1e-10):
        if abs(self.nodes[i][0]-self.nodes[j][0]-self.h)<=tol and self.nodes[i][1]==self.nodes[j][1]:
            return 'Right'
        elif abs(self.nodes[i][0]-self.nodes[j][0]+self.h)<=tol and self.nodes[i][1]==self.nodes[j][1]:
            return 'Left'
        elif self.nodes[i][0]==self.nodes[j][0] and abs(self.nodes[i][1]-self.nodes[j][1]-self.h)<=tol:
            return 'Top'
        elif self.nodes[i][0]==self.nodes[j][0] and abs(self.nodes[i][1]-self.nodes[j][1]+self.h)<=tol:
            return 'Low'
        elif self.nodes[i][0]==self.nodes[j][0] and self.nodes[i][1]==self.nodes[j][1]:
            return 'Same'
        else:
            return 'None'

class Poisson(Elliptic):


    def __init__(self, nodal_value=12, length=100,
                f=lambda x,y: 
                (2*np.pi**2/(100**2)) * np.sin(x*np.pi/100) * np.sin(y*np.pi/100)
                ):
        super().__init__(nodal_value, length)
        self.integrand_bilinear_form = lambda x, y, i, j:(
            self.phi_x(x,y,i) * self.phi_x(x,y,j) 
            + self.phi_y(x,y,i) * self.phi_y(x,y,j))
        self.f= f
        self.C=2*np.pi**2/(100**2)
    
    def _l(self,i):
        integrand=lambda x: self.phi(x[0],x[1],i)*self.f(x[0],x[1])
        h_shift=np.array([self.h,0])
        v_shift=np.array([0,self.h])
        I=0
        I+=quad(self.data46,integrand,# 1
            self.nodes[i],self.nodes[i]+h_shift,self.nodes[i]+v_shift,self.h)
        I+=quad(self.data46,integrand,# 2
            self.nodes[i],self.nodes[i]+v_shift,self.nodes[i]+v_shift-h_shift,self.h)
        I+=quad(self.data46,integrand,# 3
            self.nodes[i],self.nodes[i]+v_shift-h_shift,self.nodes[i]-h_shift,self.h)
        I+=quad(self.data46,integrand,# 4
            self.nodes[i],self.nodes[i]-h_shift,self.nodes[i]-v_shift,self.h)
        I+=quad(self.data46,integrand,# 5
            self.nodes[i],self.nodes[i]-v_shift,self.nodes[i]-v_shift+h_shift,self.h)
        I+=quad(self.data46,integrand,# 6
            self.nodes[i],self.nodes[i]-v_shift+h_shift,self.nodes[i]+h_shift,self.h)
        return I

    def _b(self):
        """Returns and store as an attribute the b vector of the final linear system."""
        self.b=np.zeros(len(self.nodes))
        for i in trange(len(self.nodes)):
            self.b[i]=self._l(i)
        return "Done!"

    def _a(self,i,j):
        integrand=lambda x: self.integrand_bilinear_form(x[0],x[1],i,j)
        h_shift=np.array([self.h,0])
        v_shift=np.array([0,self.h])
        if self.isonthe(i,j)=='Right':
            I=0
            I+=quad(self.data16,integrand,
                self.nodes[i],self.nodes[j],self.nodes[j]+v_shift,self.h)
            I+=quad(self.data16,integrand,
                self.nodes[i],self.nodes[j],self.nodes[i]-v_shift,self.h)
            return I
        elif self.isonthe(i,j)=='Left':
            I=0 
            I+=quad(self.data16,integrand,
                self.nodes[i],self.nodes[j],self.nodes[i]+v_shift,self.h)
            I+=quad(self.data16,integrand,
                self.nodes[i],self.nodes[j],self.nodes[j]-v_shift,self.h)
            return I
        elif self.isonthe(i,j)== 'Top':
            I=0 
            I+=quad(self.data16,integrand,
                self.nodes[i],self.nodes[j],self.nodes[i]-h_shift,self.h)
            I+=quad(self.data16,integrand,
                self.nodes[i],self.nodes[j],self.nodes[j]+h_shift,self.h)
            return I
        elif self.isonthe(i,j)=='Low':
            I=0 
            I+=quad(self.data16,integrand,
                self.nodes[i],self.nodes[j],self.nodes[j]-h_shift,self.h)
            I+=quad(self.data16,integrand,
                self.nodes[i],self.nodes[j],self.nodes[i]+h_shift,self.h)
            return I
        elif self.isonthe(i,j)=='Same':
            I=0
            I+=quad(self.data16,integrand,# 1
                self.nodes[i],self.nodes[i]+h_shift,self.nodes[i]+v_shift,self.h)
            I+=quad(self.data16,integrand,# 2
                self.nodes[i],self.nodes[i]+v_shift,self.nodes[i]+v_shift-h_shift,self.h)
            I+=quad(self.data16,integrand,# 3
                self.nodes[i],self.nodes[i]+v_shift-h_shift,self.nodes[i]-h_shift,self.h)
            I+=quad(self.data16,integrand,# 4
                self.nodes[i],self.nodes[i]-h_shift,self.nodes[i]-v_shift,self.h)
            I+=quad(self.data16,integrand,# 5
                self.nodes[i],self.nodes[i]-v_shift,self.nodes[i]-v_shift+h_shift,self.h)
            I+=quad(self.data16,integrand,# 6
                self.nodes[i],self.nodes[i]-v_shift+h_shift,self.nodes[i]+h_shift,self.h)
            return I
        else:
            return 0

    def _A(self):
        """Returns and store as an attribute the A matrix of the final linear system."""
        self.A=np.zeros((len(self.nodes),len(self.nodes)))
        for i in trange(len(self.nodes)): # the code uses trange to get loading bar
            for j in range(len(self.nodes)):
                self.A[i,j]=self._a(i,j)
        return "Done!"

    def _U(self):
        """Returns and store as an attribute the U vector of coefficients."""
        if isinstance(self.A,str) or isinstance(self.b,str):
            self._A()
            self._b()
        self.U=linalg.solve(self.A,self.b)
        return 'Done!'
    
    def uh(self,x,y): 
        """Approximated solution u(x,y)."""
        if isinstance(self.U,str):
            self._U()
        res=0
        for i in range(len(self.nodes)):
            res+=self.U[i]*self.phi(x,y,i)
        return res

    def u(self,x,y):
        return self.f(x,y)/self.C        
        
class Helmotz(Poisson):
    

    def __init__(self,c=lambda x,y: 4, nodal_value=12, length=100,
        dataframe_filename = 16,
        f=lambda x,y: (2*np.pi**2/100**2+1)*np.sin(x*np.pi/100)*np.sin(y*np.pi/100)
        ):
        super().__init__(nodal_value, length, dataframe_filename,f)
        self.c = c
        self.integrand_bilinear_form = lambda x,y,i,j: (
        self.phi_x(x, y, i) * self.phi_x(x, y, j) 
        + self.phi_y(x, y, i) * self.phi_y(x, y, j) 
        + self.phi(x, y, i) * self.c(x, y) )

