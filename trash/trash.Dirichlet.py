# def nodal_basis_codomain(mesh,nodal_point,h):
#     xd=mesh[0]
#     yd=mesh[1]
#     res=np.zeros(dim(xd))
#     for i in range(dim(xd)[0]):
#         for j in range(dim(xd)[1]):  
#             res[i,j]= nodal_basis(xd[i,j],yd[i,j],nodal_point,h)
#     return res