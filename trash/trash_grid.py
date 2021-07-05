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
# def _lk(self, k, f = lambda x,y: np.sin(x+y), GL_degree=100):
    #     """Linear form of the variational formulation of 2D-Poisson PDE."""
    #     integrand_linear_form=lambda x: self.phi(x[0],x[1],k)*f(x[0],x[1])
    #     return GaussLegendre1(integrand_linear_form,self.a,self.b,self.c,self.d,GL_degree)

 # def _a(self,i,j,GL_degree=100):
    #     """Bilinear form of the variational formulation of 
    #     2D-Poisson PDE with Dirichelt BCs."""
    #     integrand=lambda x,y: self.integrand_bilinear_form(x,y,i,j)
    #     return GaussLegendre1(
    #         integrand,self.support(i,j),GL_degree)
# def _A(self,GL_degree=30):
    #     """Returns and store as an attribute the A matrix of the final linear system."""
    #     self.A=np.zeros((len(self.nodes),len(self.nodes)))
    #     for i in trange(len(self.nodes)): # the code uses trange to get loading bar
    #         for j in range(len(self.nodes)):
    #             self.A[i,j]=self._a(i,j,GL_degree)
    #     return self.A

# def _l(self, k, f = lambda x,y: np.sin(x+y), GL_degree=100):
    #     """Linear form of the variational formulation of 2D-Poisson PDE."""
    #     integrand_linear_form=lambda x,y: self.phi(x,y,k)*f(x,y)
    #     return GaussLegendre1(integrand_linear_form,self.a,self.b,self.c,self.d,GL_degree)
    
# def _L(self, f = lambda x,y: np.sin(x+y), GL_degree = 100):
    #     """Returns and store as an attribute the b vector of the final linear system."""
    #     L=np.zeros(len(self.nodes))
    #     for i in range(len(self.nodes)):
    #         L[i]=self._l(i,f,GL_degree)
    #     self.L=L
    #     return L


# class Elliptic:
#     """Initialise data for solving an Elliptic PDE in a 2D square domain.

#     Contains information about the domain and the cells that will
#     be used in the finite element approximation of the domain.

#     Parameteres
#     ----------
#     nodal_value : integer
#         The desired number of nodes in the first row of nodes;
#     length : float
#         The desired lenght of the side of square domain.
#     origin : 2-sized list
#         The desired point from where the square-domain will start,


#     Attributes
#     ----------
#     nodes : np.ndarray
#         Cartesian coordinates of the interior nodes.
#     h : float
#         Diameter
#     length : float
#         Length of the side of the square-domain
#     phi : function(x : float , y : float, k : integer) -> float
#         Piecewise linear function of the finite element method
#     phi_x : function(x : float , y : float, k : integer) -> float
#         Partial derivative w.r.t x of phi
#     phi_y : function(x : float , y : float, k : integer) -> float
#         Partial derivative w.r.t y of phi
#     data16 : pd.core.frame.DataFrame
#         Barycentric coordinates and weights for 16 points triangular quadrature
#     data46 : pd.core.frame.DataFrame
#         Barycentric coordinates and weights for 46 points trinagular quadrature
#     domain : list( float, float, float, float )
#         Contains information descibing the location of the square-domain.

#     """

#     def __init__(self, nodal_value=12, length=100, origin=[0, 0]):
#         xd, yd = np.meshgrid(
#             np.linspace(origin[0], length+origin[0], nodal_value),
#             np.linspace(origin[1], length+origin[1], nodal_value))
#         self.nodes = np.transpose(
#             (np.concatenate(xd[1:-1, 1:-1]), np.concatenate(yd[1:-1, 1:-1])))
#         self.h = xd[0, 1]-xd[0, 0]
#         self.length = length
#         self.phi = lambda x, y, k: nodal_basis(x, y, self.nodes[k], self.h)
#         self.phi_x = lambda x, y, k: nodal_basis_x(x, y, self.nodes[k], self.h)
#         self.phi_y = lambda x, y, k: nodal_basis_y(x, y, self.nodes[k], self.h)
#         self.data16 = pd.read_csv("data/16.csv")
#         self.data46 = pd.read_csv("data/46.csv")
#         self.domain = [origin[0], length+origin[0],
#                        origin[1], length+origin[1]]

#     def __repr__(self):
#         return self.__class__.__name__ + " with x in " + repr(
#             self.domain[0]) + " , " + " y in "+repr(
#             self.domain[1]) + " and " + repr(len(self.nodes)) + " cells "

#     def isonthe(self, i, j, tol=1e-10):
#         """Returns a string saying where the i-th interior node is with
#         respect to the j-th one."""

#         if not -len(self.nodes) < i < len(self.nodes):
#             raise ValueError(f"{i} is not between {-len(self.nodes)}"
#                              f"and {len(self.nodes)}")
#         if not -len(self.nodes) < j < len(self.nodes):
#             raise ValueError(f"{j} is not between {-len(self.nodes)}"
#                              f" and {len(self.nodes)}")
#         if self.nodes[i][1] == self.nodes[j][1]:
#             if abs(self.nodes[i][0]-self.nodes[j][0] - self.h) <= tol:
#                 return 'Right'
#             elif abs(self.nodes[i][0]-self.nodes[j][0]+self.h) <= tol:
#                 return 'Left'
#             elif self.nodes[i][0] == self.nodes[j][0]:
#                 return 'Same'
#             else:
#                 return 'None'
#         elif self.nodes[i][0] == self.nodes[j][0]:
#             if abs(self.nodes[i][1]-self.nodes[j][1] - self.h) <= tol:
#                 return 'Top'
#             elif abs(self.nodes[i][1]-self.nodes[j][1]+self.h) <= tol:
#                 return 'Low'
#             else:
#                 return 'None'
#         elif abs(self.nodes[i][1] - self.nodes[j][1] - self.h) <= tol:
#             if abs(self.nodes[i][0]-self.nodes[j][0] + self.h) <= tol:
#                 return 'TopLeft'
#             else:
#                 return 'None'
#         elif abs(self.nodes[i][1] - self.nodes[j][1] + self.h) <= tol:
#             if abs(self.nodes[i][0]-self.nodes[j][0] - self.h) <= tol:
#                 return 'LowRight'
#             else:
#                 return 'None'
#         else:
#             return 'None'


# class Poisson(Elliptic):
#     """Subclass of Elliptic for solving Poisson equation.

#         Contains methods for solving Poisson equation on 2D Square Domain
#         by finite element method with piecewise linear functions.

#         Parameteres
#         ----------
#         nodal_value : integer
#             The desired number of nodes in the first row of nodes;
#         length : float
#             The desired lenght of the side of square domain.
#         origin : 2-sized list
#             The desired point from where the square-domain will start.
#         f: function(x : float, y: float) -> float
#             The f(x,y) function of the Poisson equation.
#         u: function(x : float, y: float) -> float (Optional)
#             The u(x,y) analytical solution of the Poisson equation

#         Attributes
#         ----------
#         nodes : np.ndarray
#             Cartesian coordinates of the interior nodes.
#         h : float
#             Diameter
#         length : float
#             Length of the side of the square-domain
#         phi : function(x : float , y : float, k : integer) -> float
#             Piecewise linear function of the finite element method
#         phi_x : function(x : float , y : float, k : integer) -> float
#             Partial derivative w.r.t x of phi
#         phi_y : function(x : float , y : float, k : integer) -> float
#             Partial derivative w.r.t y of phi
#         data16 : pd.core.frame.DataFrame
#             Barycentric coordinates and weights for
#             16 points triangular quadrature.
#         data46 : pd.core.frame.DataFrame
#             Barycentric coordinates and weights for
#             46 points trinagular quadrature.
#         domain : list( float, float, float, float )
#             Contains information descibing the location of the square-domain.
#         A : str or np.ndarray
#             Stiffness matrix
#         b : str or np.ndarray
#             Right hand side of linear system with stiffness matrix.
#         U : str or np.ndarray
#             Solution array of linear system.
#         L2error : float (Optional)
#             L2norm of the error between analytical
#             solution and Galerkin approximation.
#         integrand_bilinear_form: function(x : float, y : float ,
#                                           i : int, j : int) -> float
#             The integrand of the bilenar form at nodes i and j.
#         f = function (x : float, y : float) -> float
#             The f(x,y) function of the Poisson equation.
#         u = function (x : float, y : float) -> float (Optional)
#             The u(x,y) analytical solution of the Poisson equation.
#         """

#     def __init__(self, nodal_value=12, length=100, origin=[0, 0],
#                  f=lambda x, y:
#                  (2*np.pi**2/(100**2)) * np.sin(x*np.pi/100)
#                  * np.sin(y*np.pi/100),
#                  u=lambda x, y:
#                  np.sin(x*np.pi/100) * np.sin(y*np.pi/100)
#                  ):
#         super().__init__(nodal_value, length, origin)
#         self.A = "Need to run the _A() method."
#         self.b = "Need to run the _L() method."
#         self.U = "Need to run the _U() method."
#         self.L2error = "Need to run the _L2error() method."
#         self.integrand_bilinear_form = lambda x, y, i, j: (
#             self.phi_x(x, y, i) * self.phi_x(x, y, j)
#             + self.phi_y(x, y, i) * self.phi_y(x, y, j))
#         self.f = f
#         self.u = u

#     def __repr__(self):
#         return self.__class__.__name__ + " with x in " + repr(
#             self.domain[0]) + " , " + "y in "+repr(
#             self.domain[1]) + " , " + repr(
#             len(self.nodes)) + " cells " + " and f = " + repr(
#             self.f)

#     def _l(self, i):
#         """Approximate the linear form at interior node 'i' """
#         def integrand(x):
#             return self.phi(x[0], x[1], i)*self.f(x[0], x[1])
#         h_shift = np.array([self.h, 0])
#         v_shift = np.array([0, self.h])
#         integral = 0
#         integral += quad(self.data46, integrand,  # 1
#                          self.nodes[i], self.nodes[i]+h_shift,
#                          self.nodes[i]+v_shift, self.h)
#         integral += quad(self.data46, integrand,  # 2
#                          self.nodes[i], self.nodes[i]+v_shift,
#                          self.nodes[i]+v_shift-h_shift, self.h)
#         integral += quad(self.data46, integrand,  # 3
#                          self.nodes[i], self.nodes[i]+v_shift-h_shift,
#                          self.nodes[i]-h_shift, self.h)
#         integral += quad(self.data46, integrand,  # 4
#                          self.nodes[i], self.nodes[i]-h_shift,
#                          self.nodes[i]-v_shift, self.h)
#         integral += quad(self.data46, integrand,  # 5
#                          self.nodes[i], self.nodes[i]-v_shift,
#                          self.nodes[i]-v_shift+h_shift, self.h)
#         integral += quad(self.data46, integrand,  # 6
#                          self.nodes[i], self.nodes[i]-v_shift+h_shift,
#                          self.nodes[i]+h_shift, self.h)
#         return integral

#     def _b(self):
#         """Returns and store as an attribute the b vector."""
#         self.b = np.zeros(len(self.nodes))
#         for i in trange(len(self.nodes)):
#             self.b[i] = self._l(i)
#         return "Done!"

#     def _a(self, i, j):
#         """Approximate the bilinear form only on the correspoding support
#            of the coupled nodes (i,j)."""
#         def integrand(x):
#             return self.integrand_bilinear_form(x[0], x[1], i, j)
#         h_shift = np.array([self.h, 0])
#         v_shift = np.array([0, self.h])
#         if self.isonthe(i, j) == 'Right':
#             integral = 0
#             integral += quad(self.data16, integrand,
#                              self.nodes[i], self.nodes[j],
#                              self.nodes[j]+v_shift, self.h)
#             integral += quad(self.data16, integrand,
#                              self.nodes[i], self.nodes[j],
#                              self.nodes[i]-v_shift, self.h)
#             return integral
#         elif self.isonthe(i, j) == 'Left':
#             integral = 0
#             integral += quad(self.data16, integrand,
#                              self.nodes[i], self.nodes[j],
#                              self.nodes[i]+v_shift, self.h)
#             integral += quad(self.data16, integrand,
#                              self.nodes[i], self.nodes[j],
#                              self.nodes[j]-v_shift, self.h)
#             return integral
#         elif self.isonthe(i, j) == 'Top':
#             integral = 0
#             integral += quad(self.data16, integrand,
#                              self.nodes[i], self.nodes[j],
#                              self.nodes[i]-h_shift, self.h)
#             integral += quad(self.data16, integrand,
#                              self.nodes[i], self.nodes[j],
#                              self.nodes[j]+h_shift, self.h)
#             return integral
#         elif self.isonthe(i, j) == 'Low':
#             integral = 0
#             integral += quad(self.data16, integrand,
#                              self.nodes[i], self.nodes[j],
#                              self.nodes[j]-h_shift, self.h)
#             integral += quad(self.data16, integrand,
#                              self.nodes[i], self.nodes[j],
#                              self.nodes[i]+h_shift, self.h)
#             return integral
#         elif self.isonthe(i, j) == 'Same':
#             integral = 0
#             integral += quad(self.data16, integrand,  # 1
#                              self.nodes[i], self.nodes[i]+h_shift,
#                              self.nodes[i]+v_shift, self.h)
#             integral += quad(self.data16, integrand,  # 2
#                              self.nodes[i], self.nodes[i]+v_shift,
#                              self.nodes[i]+v_shift-h_shift, self.h)
#             integral += quad(self.data16, integrand,  # 3
#                              self.nodes[i], self.nodes[i]+v_shift-h_shift,
#                              self.nodes[i]-h_shift, self.h)
#             integral += quad(self.data16, integrand,  # 4
#                              self.nodes[i], self.nodes[i]-h_shift,
#                              self.nodes[i]-v_shift, self.h)
#             integral += quad(self.data16, integrand,  # 5
#                              self.nodes[i], self.nodes[i]-v_shift,
#                              self.nodes[i]-v_shift+h_shift, self.h)
#             integral += quad(self.data16, integrand,  # 6
#                              self.nodes[i], self.nodes[i]-v_shift+h_shift,
#                              self.nodes[i]+h_shift, self.h)
#             return integral
#         else:
#             return 0

#     def _A(self):
#         """Returns and store as an attribute the A matrix. """
#         self.A = np.zeros((len(self.nodes), len(self.nodes)))
#         for i in trange(len(self.nodes)):  # trange to get loading bar
#             for j in range(len(self.nodes)):
#                 self.A[i, j] = self._a(i, j)
#         return "Done!"

#     def _U(self):
#         """Store as an attribute the U vector of coefficients."""
#         if isinstance(self.A, str) or isinstance(self.b, str):
#             self.A = np.zeros((len(self.nodes), len(self.nodes)))
#             self.b = np.zeros(len(self.nodes))
#             for i in trange(len(self.nodes)):  # trange to get loading bar
#                 self.b[i] = self._l(i)
#                 for j in range(len(self.nodes)):
#                     self.A[i, j] = self._a(i, j)
#         self.U = linalg.solve(self.A, self.b)
#         return 'Done!'

#     def uh(self, x, y):
#         """Approximated solution evaluated at cartesian coordinates (x,y)."""
#         if isinstance(self.U, str):
#             self._U()
#         res = 0
#         for i in range(len(self.nodes)):
#             res += self.U[i]*self.phi(x, y, i)
#         return res

#     def error(self, GLdegree=7, norm='L2norm'):
#         if norm == 'L2norm':
#             if isinstance(self.L2error, str):
#                 self.uh(self.length/2, self.length/2)

#             def u_h(x, y):
#                 return self.uh(x, y)

#             def u_(x, y):
#                 return self.u(x, y)
#             self.L2error = L2error(u_h, u_, self.domain[0],
#                                    self.domain[1], self.domain[2],
#                                    self.domain[3], GLdegree)
#             return self.L2error