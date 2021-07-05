
def nodal_basis(x, y, nodal_point, h):
    """Piecewise linear basis function."""
    x_i = nodal_point[0]
    y_j = nodal_point[1]
    if (y-y_j) >= 0 and (x-x_i) >= 0 and (y-y_j) <= -(x-x_i)+h:
        return 1 - (x-x_i)/h - (y-y_j)/h  # 1
    elif -h <= (x-x_i) <= 0 and (y-y_j) >= -(x-x_i) and h >= (y-y_j) >= 0:
        return 1 - (y-y_j)/h  # 2
    elif -h <= (x-x_i) <= 0 and (y-y_j) <= -(x-x_i) and h >= (y-y_j) >= 0:
        return 1 - (x_i - x)/h  # 3
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h:
        return 1 - (x_i-x)/h - (y_j-y)/h  # 4
    elif (h+x_i) >= x >= x_i and (y-y_j) <= -(x-x_i) and (-h+y_j) <= y <= y_j:
        return 1 - (y_j - y)/h  # 5
    elif (h+x_i) >= x >= x_i and (y-y_j) >= -(x-x_i) and (-h+y_j) <= y <= y_j:
        return 1 - (x-x_i)/h  # 6
    elif x == x_i and y == y_j:
        return 1
    else:
        return 0


def nodal_basis_y(x, y, nodal_point, h):
    """Partial derivative with respect of x of
       piecewise linear basis function. """
    x_i = nodal_point[0]
    y_j = nodal_point[1]
    if (y-y_j) >= 0 and (x-x_i) >= 0 and (y-y_j) <= -(x-x_i)+h:
        return -1/h  # 1
    elif -h <= (x-x_i) <= 0 and (y-y_j) >= -(x-x_i) and h >= (y-y_j) >= 0:
        return -1/h  # 2
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h:
        return 1/h  # 4
    elif (h+x_i) >= x >= x_i and (y-y_j) <= -(x-x_i) and (-h+y_j) <= y <= y_j:
        return 1/h  # 5
    else:
        return 0


def nodal_basis_x(x, y, nodal_point, h):
    """Partial derivative with respect of y of
       piecewise linear basis function. """
    x_i = nodal_point[0]
    y_j = nodal_point[1]
    if (y-y_j) >= 0 and (x-x_i) >= 0 and (y-y_j) <= -(x-x_i)+h:
        return -1/h   # 1
    elif -h <= (x-x_i) <= 0 and (y-y_j) <= -(x-x_i) and h >= (y-y_j) >= 0:
        return 1/h  # 3
    elif (x-x_i) <= 0 and (y-y_j) <= 0 and (y-y_j) >= -(x-x_i)-h:
        return 1/h  # 4
    elif (h+x_i) >= x >= x_i and (y-y_j) >= -(x-x_i) and (-h+y_j) <= y <= y_j:
        return -1/h  # 6
    else:
        return 0