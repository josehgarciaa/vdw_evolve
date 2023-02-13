"""

"""

import numpy as np
from numpy.linalg import det, inv

from .utils import t_cel1t_cel2, strained_proces


def base_fit_function(params, cel1, cel2, strain_boundary):
    """
    This creates a cost that the annealing process will minimize.
    The cost is proportional to the area of the supercell
    and the length of the new based vectors.
    :param params: tA parameters
    :param cel1:
    :param cel2:
    :param strain_boundary:
    :return:
    """
    t_cel1, strain_t_cel2 = t_cel1t_cel2(params, cel1, cel2)
    super_c_1 = np.dot(t_cel1, cel1)

    # It's a chance that a solution with t_cel2 integer does not exist. So if the t_cel2 has float values,
    # we will decompose it in a new t_cel2 integer and a strain(diagonal matrix), which multiply and reconstruct
    # the initial matrix.
    strain = strained_proces(strain_t_cel2, strain_boundary)
    t_cel2 = np.dot(inv(strain), strain_t_cel2)

    # area of the new cel area = det(cel1)*det(t_cel1)
    super_cel_area = det(super_c_1) * det(super_c_1)  # minimum but bigger than 0

    # t_cel1 integer # cost for t_cel1 not be an integer
    cons = 99999
    round_cost = 0
    for row in t_cel2:
        for e in row:
            z_ero1 = (round(e) - e) * (round(e) - e)
            round_cost += z_ero1
    round_cost = round_cost * cons

    # cost for minimizing the new vectors length
    # ad 1  0,1**n <1
    length_cost = (1+(t_cel1[0][0] ** 2 + t_cel1[0][1] ** 2) + (t_cel1[1][0] ** 2 + t_cel1[1][1] ** 2)) * 100

    f = partial_cost_shape(super_cel_area) + round_cost ** 3 + length_cost * 2

    return f


def partial_cost_shape(x, up=99999999):
    """
    This function will take the area of the new cell as a parameter
    it is built in such that it will take the minimum value when
    the area is close to 1 and will go to inf in 0 and infinity.

    :param x: float  super_cel_area
    :param up: float  modul value of the peek in 1.
    :return: float
    """
    tr_x = (x + 0.5)

    if x < 1:
        res = up * (1 / (1 + np.sqrt((x ** 2)))) * np.sin(tr_x * np.pi) + 1 / (x + 0.000000000000001)
    else:
        res = (x - 1) ** 2 - up * (1 / (1 + x ** 2)) + 1 / (x + 0.00000000000001)
    return res




###



def rectangle_fit_function(params, cel1, cel2, strain_boundary):
    """
    This creates a cost that the annealing process will minimize.
    The cost is proportional to the area of the supercell
    and the length of the new based vectors.
    :param params: tA parameters
    :param cel1:
    :param cel2:
    :param strain_boundary:
    :return:
    """
    t_cel1, strain_t_cel2 = t_cel1t_cel2(params, cel1, cel2)
    super_c_1 = np.dot(t_cel1, cel1)



    # cost for minimizing the new vectors length
    # thry with this constant 
    length_cost = (1+(t_cel1[0][0] ** 2 + t_cel1[0][1] ** 2) + (t_cel1[1][0] ** 2 + t_cel1[1][1] ** 2))

    f =  (1+super_c_1[0][1]**2 +super_c_1[1][0]**2)**5*100 + length_cost * 2
    # print("da")

    # print(f)
    return f
