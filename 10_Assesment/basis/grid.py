import numpy as np 
import numbers
from basis.basis import basis

class SqGrid:

    def __init__(self,l=100,n=4,m=25):
        
        if not isinstance(n,numbers.Integral):
            raise TypeError(f"{n} is not an integer")
        
        xd , yd = np.meshgrid(np.linspace(0,l,n),np.linspace(0,l,n))
        hypotenuse = np.sqrt(xd[0,1]**2 + yd[1,0]**2)
        nodes = np.transpose((np.concatenate(xd[1:-1,1:-1]),np.concatenate(yd[1:-1,1:-1])))
        #notice: Number of nodes = (n-2)(n-2)
        xv,yv = np.meshgrid(np.linspace(0,l,round(n*m)),np.linspace(0,l,round(n*m)))
        

        self.dim = (len(xv),len(yv[0]))
        self.h = np.sin(np.pi/4)*hypotenuse
        self.nodes = nodes
        self.mesh = xv,yv

    def evaluate(self,k):
        res=basis(self.mesh,self.nodes[k],self.h)
        return res

    def check(self,k):
        res = self.evaluate(k)
        if 0.965<=res[round(self.nodes[k][1]),round(self.nodes[k][0])]<=1:
            return f"({round(self.nodes[k][0])},{round(self.nodes[k][1])}) seems correct."
        elif 0.94<=res[round(self.nodes[k][1]),round(self.nodes[k][0])]<0.965:
            raise Warning(f"({round(self.nodes[k][0])},{round(self.nodes[k][1])}) might not be correct.")
        else:
            raise EstimationError("Try some coccobello.")

    def clean_eva(self,k):
        res=self.evaluate(k)
        for i in range(len(res)):
            for j in range(len(res[0])):
                if res[i,j]<0:
                    res[i,j]=0
        return res


    def machine_errors(self,k):
        res=self.evaluate(k)
        logic=[]
        for i in range(len(res)):
            for j in range(len(res[0])):
                if res[i,j]<0:
                    logic+=f"{(i,j)}"
        return logic

class EstimationError(Exception):
    pass

    

    




  
    

