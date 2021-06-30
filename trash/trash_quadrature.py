#def GaussLegendre3(phi,f,ksq,a,b,c,d,n=7):
#     [xi,w]=p_roots(n+1)
#     x=1/2*(1-xi)*a+1/2*(1+xi)*b
#     y=1/2*(1-xi)*c+1/2*(1+xi)*d
#     G=np.zeros(((n+1),(n+1)))
#     for i in range(n+1):
#         for j in range(n+1):
#             G[i,j]= w[i]*w[j]*phi(x[i],y[j])*f(x[i],y[j])*ksq(x[i],y[i])
#     G=np.sum(G)*(b-a)/2*(d-c)/2
#     return G

# def GaussLegendre(phi,f,a,b,c,d,n=7):
#     [xi,w]=p_roots(n+1)
#     x=1/2*(1-xi)*a+1/2*(1+xi)*b
#     y=1/2*(1-xi)*c+1/2*(1+xi)*d
#     G=np.zeros(((n+1),(n+1)))
#     for i in range(n+1):
#         for j in range(n+1):
#             G[i,j]= w[i]*w[j]*phi(x[i],y[j])*f(x[i],y[j])
#     G=np.sum(G)*(b-a)/2*(d-c)/2
#     return G

# def GaussLegendre1(function,interval,n=7):
#     [xi,w]=p_roots(n+1)
#     x=1/2*(1-xi)*interval[0]+1/2*(1+xi)*interval[1]
#     y=1/2*(1-xi)*interval[2]+1/2*(1+xi)*interval[3]
#     G=np.zeros(((n+1),(n+1)))
#     for i in range(n+1):
#         for j in range(n+1):
#             G[i,j]= w[i]*w[j]*function(x[i],y[j])
#     return np.sum(G)*(interval[1]-interval[0])/2*(interval[3]-interval[2])/2  