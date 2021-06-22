from basis.basis import finite_basis_grid
import numpy as np

xd,yd=np.meshgrid(np.linspace(0,1,100),np.linspace(0,1,100))

res11=finite_basis_grid(xd,yd,0.2,0.2,np.sqrt(0.2**2+0.2**2))
