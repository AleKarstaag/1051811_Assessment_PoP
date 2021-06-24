
def nodal_basis(x,y,nodal_point,h):
    """Finite basis function."""
    x_i=nodal_point[0]
    y_j=nodal_point[1]
    if (y-y_j)>=0 and (x-x_i)>=0 and (y-y_j) <= -(x-x_i)+h: #1
        return 1 - (x-x_i)/h - (y-y_j)/h
    elif -h <= (x-x_i) <= 0 and (y-y_j) >= -(x-x_i) and h>= (y-y_j) >= 0: #2
        return 1 - (y-y_j)/h
    elif -h <= (x-x_i) <= 0 and (y-y_j) <= -(x-x_i) and h>=(y-y_j)>=0: #3
        return 1 - (x_i - x)/h
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h : #4  
        return 1 - (x_i-x)/h - (y_j-y)/h
    elif (h+x_i)>=x>= x_i and (y-y_j) <= -(x-x_i) and (-h+y_j)<=y<=y_j: #5
        return 1 - (y_j - y)/h
    elif (h+x_i)>=x>=x_i and (y-y_j) >= -(x-x_i) and (-h+y_j) <=y<=y_j: #6
        return 1 - (x-x_i)/h
    elif x==x_i and y==y_j:
        return 1    
    else:
        return 0

def nodal_basis_y(x,y,nodal_point,h):
    """Finite basis function."""
    x_i=nodal_point[0]
    y_j=nodal_point[1]
    if (y-y_j)>=0 and (x-x_i)>=0 and (y-y_j) <= -(x-x_i)+h: #1
        return -1/h
    elif -h <= (x-x_i) <= 0 and (y-y_j) >= -(x-x_i) and h>= (y-y_j) >= 0: #2
        return -1/h 
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h : #4  
        return 1/h
    elif (h+x_i)>=x>= x_i and (y-y_j) <= -(x-x_i) and (-h+y_j)<=y<=y_j: #5
        return 1/h
    else:
        return 0

def nodal_basis_x(x,y,nodal_point,h): 
    """Finite basis function."""
    x_i=nodal_point[0]
    y_j=nodal_point[1]
    if (y-y_j)>=0 and (x-x_i)>=0 and (y-y_j) <= -(x-x_i)+h: #1
        return -1/h
    elif -h <= (x-x_i) <= 0 and (y-y_j) <= -(x-x_i) and h>=(y-y_j)>=0: #3
        return 1/h
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h : #4  
        return 1/h
    elif (h+x_i)>=x>=x_i and (y-y_j) >= -(x-x_i) and (-h+y_j) <=y<=y_j: #6
        return -1/h
    else:
        return 0

