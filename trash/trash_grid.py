#   class Elliptic:
    #def _evaluate(self,k):
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
    
# class Class1:
#     def __init__(self,value1,value2):
#         self.attr1=value1
#         self.attr2=value1*2
#         self.attr3=value1+value2

# class Class2(Class1):
#     def __init__(self, value1,value2,extra_val1):
#         super().__init__(value1,value2)
#         self.extra=extra_val1

# phi_i_x=lambda x,y: nodal_basis_x(x,y,self.nodes[i],self.h)
# phi_j_x=lambda x,y: nodal_basis_x(x,y,self.nodes[j],self.h)
# phi_i_y=lambda x,y: nodal_basis_y(x,y,self.nodes[i],self.h)
# phi_j_y=lambda x,y: nodal_basis_y(x,y,self.nodes[j],self.h)

        #     xv, yv = np.meshgrid(
        #     np.linspace(0, length, round(n*m)),
        #     np.linspace(0, length, round(n*m)))